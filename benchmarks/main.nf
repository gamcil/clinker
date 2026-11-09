nextflow.enable.dsl = 2

params.fixtures = "${projectDir}/fixtures/*.tar.gz"
params.repeats = 3
params.threads = 4
params.identity = 0.3
params.kmer_size = 5
params.min_shared_kmers = 8
params.container = 'ghcr.io/gamcil/clinker-benchmarks:0.1.0'

process EXTRACT_FIXTURE {
    tag { id }
    container params.container
    maxForks 1
    cache 'deep'

    input:
    tuple val(id), val(role), val(clusters), val(proteins), val(sha256), path(archive)

    output:
    tuple val(id), val(role), val(clusters), val(proteins), path(id)

    script:
    """
    echo '${sha256}  ${archive}' | sha256sum --check --status
    mkdir -p '${id}'
    tar -xzf '${archive}' --strip-components=1 -C '${id}'
    """
}

process BENCHMARK {
    tag { "${id} ${mode} repeat ${repeat}" }
    container params.container
    cpus params.threads
    maxForks 1
    publishDir params.outdir, mode: 'copy'

    input:
    tuple val(id), val(role), val(clusters), val(proteins), path(fixture), val(mode), val(repeat)

    output:
    path("${id}.${mode}.${repeat}.json"), emit: records
    path("${id}.${mode}.${repeat}.log"), emit: logs
    path("${id}.${mode}.${repeat}.links"), emit: links

    script:
    """
    python /opt/benchmark/measure.py \
      --fixture '${fixture}' \
      --role '${role}' \
      --clusters ${clusters} \
      --proteins ${proteins} \
      --mode '${mode}' \
      --repeat ${repeat} \
      --threads ${task.cpus} \
      --identity ${params.identity} \
      --kmer-size ${params.kmer_size} \
      --min-shared-kmers ${params.min_shared_kmers} \
      --container '${params.container}' \
      --output '${id}.${mode}.${repeat}'
    """
}

process SUMMARIZE {
    container params.container
    publishDir params.outdir, mode: 'copy'

    input:
    path(records)

    output:
    path('summary.csv')

    script:
    """
    python /opt/benchmark/summarize.py ${records.join(' ')} summary.csv
    """
}

workflow {
    fixture_metadata = file("${projectDir}/fixtures.tsv")
        .readLines()
        .drop(1)
        .collectEntries { line ->
            def row = line.split('\\t')
            [(row[0]): [role: row[1], clusters: row[2] as Integer,
                        proteins: row[3] as Integer, sha256: row[4]]]
        }

    archives = Channel
        .fromPath(params.fixtures, checkIfExists: true)
        .map { archive ->
            def id = archive.name.replaceFirst(/\.tar\.gz$/, '')
            def metadata = fixture_metadata[id]
            if (metadata == null) error("No fixtures.tsv entry for ${id}")
            tuple(id, metadata.role, metadata.clusters, metadata.proteins,
                  metadata.sha256, archive)
        }

    extracted_fixtures = EXTRACT_FIXTURE(archives).collect()

    jobs = extracted_fixtures.flatMap { fixtures ->
        fixtures.collectMany { id, role, clusters, proteins, directory ->
            ['python', 'rust', 'rust-kmer'].collectMany { mode ->
                (1..(params.repeats as int)).collect { repeat ->
                    tuple(id, role, clusters, proteins, directory, mode, repeat)
                }
            }
        }
    }

    runs = BENCHMARK(jobs)
    SUMMARIZE(runs.records.collect())
}
