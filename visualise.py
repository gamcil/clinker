#!/usr/bin/env python3

"""
Module for generating the cluster comparison figure.

Cameron Gilchrist
"""

from classes import Cluster, Gene


class Figure:
    """ SVG Figure constructor class.
    """

    # define the default colourscheme
    # taken from colorbrewer2 - paired n12
    default_colours = {
        'PKS': '#a6cee3',
        'NRPS': '#1f78b4',
        'PKS-NRPS': '#33a02c',
        'P450': '#b2df8a',
        'Oxidase': '#ffffff',
        'Regulator': '#fb9a99',
        'Hydrolase': '#e31a1c',
        'Reductase': '#fdbf6f6',
        'Phosphatase': '#ff7f00',
        'Transporter': '#cab2d6',
        'Dehalogenase': '#6a3d9a',
        'Dehydrogenase': '#ffff99',
        'Terpene synthase': '#b15926',
        'Hypothetical': '#d3d3d3'
    }

    def __init__(self,
                 colours=None,
                 vertical_spacing=None,
                 locus_spacing=None,
                 width=None,
                 aligner=None,
                 tip_length=None,
                 wing_height=None,
                 gene_height=None
                 ):
        self.aligner = aligner  # ClusterAligner object
        if not self.aligner:
            raise ValueError(
                'A Figure object must be instantiated with a ClusterAligner'
            )

        self.colours = colours if colours else self.default_colours

        # Vertical spacing between clusters
        self.cluster_spacing = vertical_spacing if vertical_spacing else 40

        # Horizontal spacing between cluster loci
        self.locus_spacing = locus_spacing if locus_spacing else 20

        # Pixel width of the figure
        self.width = width if width else 1000
        self.tip_length = tip_length if tip_length else 8
        self.wing_height = wing_height if wing_height else 0
        self.gene_height = gene_height if gene_height else 12

    def render(self, order):
        """ Render the complete SVG of this figure.
        """
        if not order:
            raise ValueError('A cluster order must be specified')
        display_order = self.aligner.compute_display_order(order)

        # Find largest cluster and configure figure scaling
        largest = self.aligner.find_largest_cluster()
        spaced_length = largest.spaced_length(self.locus_spacing)
        scale_factor = self.width / spaced_length

        attachment_points = {}
        previous_cluster = None
        cluster_maps = []
        tracks = ''  # draw these in order for correct z index
        links = ''
        genes = ''

        for count, cluster_name in enumerate(display_order):  # cluster names
            cluster = self.aligner.get_cluster(cluster_name)
            cluster_maps.append(
                cluster.map_cluster(locus_spacing=self.locus_spacing,
                                    scale_factor=scale_factor)
            )

            # Get the alignment to the previous cluster
            if count > 0:
                alignment = self.aligner.get_alignment(
                    display_order[count - 1],
                    cluster_name
                )

                first, last = None, None
                previous_total = len(previous_cluster.genes)

                for index in range(previous_total):
                    lower = previous_cluster.genes[index]
                    upper = previous_cluster.genes[previous_total - index - 1]

                    if first and last:
                        break

                    if not first and lower.name in attachment_points and \
                            alignment.find_gene(lower.name):
                        first = attachment_points[lower.name][0]

                    if not last and upper.name in attachment_points and \
                            alignment.find_gene(upper.name):
                        last = attachment_points[upper.name][2]

                x_offset = (
                    (abs(last - first) / 2) + first -
                    find_midpoint(alignment, cluster_maps[-1])
                )
            else:
                alignment, x_offset = None, 0

            for locus_name, locus in cluster_maps[-1].items():
                y_cluster = self.cluster_spacing * count

                tracks += (
                    '<line x1="{}" y1="{mid}" x2="{}" y2="{mid}" '
                    'stroke="#d3d3d3" stroke-width="2" />'
                ).format(
                    x_offset + int(locus['start']),
                    x_offset + int(locus['end']),
                    mid=y_cluster + self.gene_height / 2
                )

                for gene_name, location in locus['genes'].items():
                    start, end, strand, function = location
                    polygon, attachments = render_gene(
                        gene_name,
                        start, end, strand,
                        self.colours[function],
                        tip_length=self.tip_length,
                        wing_height=self.wing_height,
                        gene_height=self.gene_height,
                        horizontal_offset=x_offset,
                        vertical_offset=y_cluster
                    )

                    attachment_points[gene_name] = attachments
                    genes += polygon

            previous_cluster = cluster
            if not alignment:
                continue

            for link in alignment.links:
                query = attachment_points[link.query]
                target = attachment_points[link.target]

                if query[1] < target[1]:  # check vertical order
                    points = [*query, *target[2:4], *target[0:2]]
                else:
                    points = [*target, *query[2:4], *query[0:2]]

                points = [int(point) for point in points]

                #    A--------------B
                #   /              /
                #  /              /
                # D--------------C

                links += (
                    f'<polygon id="link_{link.query}_{link.target}"'
                    ' points="{},{},{},{},{},{},{},{}"'
                    ' fill="#000000" fill-opacity="{}" />'
                ).format(*points, link.identity)

        return f'<svg width="{self.width}">{tracks}{links}{genes}</svg>'


def find_midpoint(alignment, cluster_map):
    """ Find the midpoint of an alignment area in a cluster.

        i.e. |........=..===.==..................|   <-- cluster
                     |=..===.==|  <-- alignment area
    """
    first, last = alignment.links[0], alignment.links[-1]

    # Re-order if necessary; figure out if this cluster has
    # queries or targets from the alignment
    for locus in cluster_map.values():

        if first.query in locus['genes']:
            alignment_start = locus['genes'][first.query][0]
        if last.query in locus['genes']:
            alignment_end = locus['genes'][last.query][1]

        if first.target in locus['genes']:
            alignment_start = locus['genes'][first.target][0]
        if last.target in locus['genes']:
            alignment_end = locus['genes'][last.target][1]

    if not alignment_start or not alignment_end:
        return None
    return alignment_start + (alignment_end - alignment_start) / 2


def render_gene(name, start, end, strand, colour,
                tip_length=10, wing_height=0, gene_height=10,
                horizontal_offset=0, vertical_offset=0):
    """ Generate a scaled SVG representation of a gene.

        Arguments:
             gene (obj): Gene object to draw
           colour (str): HEX code to use for the fill of the
                         polygon
       tip_length (int): length of the arrow tip
      wing_height (int): height of the 'wing' of each arrow tip_length
      gene_height (int): height of the gene body

        Returns:
            polygon (str): SVG representation of the Gene
 attachment_points (list): coordinates of homology link attachment points

        A gene arrow is drawn as an SVG polygon with points A-G.

                       C           B
                       |\         /|
            A----------B \       / C--------D
            |             \     /           |
            |              D   A            |
            |             /     \           |
            G----------F /       \ F--------E
                       |/         \|
                       E           G
    """
    total_length = end - start
    arrow_body = total_length - tip_length
    y_midpoint = gene_height / 2
    y_bottom = vertical_offset + gene_height

    if strand == 1:  # forward strand genes
        tip_start = horizontal_offset + start + arrow_body

        coordinates = [  # from A to G...
            horizontal_offset + start, vertical_offset,
            tip_start, vertical_offset,
            tip_start, vertical_offset - wing_height,
            horizontal_offset + start + total_length,
            vertical_offset + y_midpoint,
            tip_start, y_bottom + wing_height,
            tip_start, y_bottom,
            horizontal_offset + start, y_bottom
        ]

    else:  # reverse strand
        tip_start = horizontal_offset + start + tip_length
        body_end = horizontal_offset + start + total_length

        coordinates = [
            horizontal_offset + start, vertical_offset + y_midpoint,
            tip_start, vertical_offset - wing_height,
            tip_start, vertical_offset,
            body_end, vertical_offset,
            body_end, y_bottom,
            tip_start, y_bottom,
            tip_start, y_bottom + wing_height,
        ]

    # Any homology links drawn will connect to these
    attachment_points = [
        horizontal_offset + start, vertical_offset + y_midpoint,
        horizontal_offset + start + total_length, vertical_offset + y_midpoint
    ]

    # Concatenate coordinates to use as 'points' arg, make polygon
    points = ','.join(str(int(c)) for c in coordinates)
    polygon = (f'<text x="{horizontal_offset + start}" y="{vertical_offset}"'
               f' font-size="12">{name}</text>\n'
               f'<polygon id="{name}" points="{points}"'
               f' fill="{colour}" stroke="black" stroke-width="1.5"/>\n')

    return polygon, attachment_points


def parse_colour_file(handle):
    """ Parse a custom colourscheme file.

        This module colours genes based on the following functions:
            Cytochrome P450
            Hydrolase
            Polyketide synthase (PKS)
            Non-ribosomal peptide synthase (NRPS)
            PKS-NRPS hybrids
            Terpene synthase
            Reductase
            Dehydrogenase
            Dehalogenase
            Phosphatase
    """
    colours = {}

    valid_functions = {
        'Cytochrome P450', 'Hydrolase', 'PKS',
        'NRPS', 'PKS-NRPS', 'Terpene synthase',
        'Reductase', 'Dehydrogenase',
        'Dehalogenase', 'Phosphatase'
    }

    for line in handle:
        try:
            function, colour = line.split(',')
        except ValueError as exc:
            # Not enough values to unpack; invalid split
            raise ValueError(
                'Custom colour file must be a 2 column CSV'
                ' e.g. P450,#FFFFFF'
            ) from exc

        colour_length = len(colour)
        if not colour.startswith('#') or colour_length != 7:
            raise ValueError(
                'Colours must be given as valid, 6 digit'
                ' HEX codes starting with #, e.g. #FFFFFF'
            )

        if function not in valid_functions:
            raise ValueError(
                'Invalid biosynthetic function supplied',
                valid_functions
            )

        colours[function] = colour

    return colours
