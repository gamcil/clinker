FROM rust:1.88.0-slim-bookworm@sha256:38bc5a86d998772d4aec2348656ed21438d20fcdce2795b56ca434cf21430d89 AS rust-build

ARG CLINKER_CLI_VERSION=0.1.0
RUN cargo install clinker-cli \
    --version "${CLINKER_CLI_VERSION}" \
    --locked \
    --root /opt/clinker-cli

FROM node:22-bookworm AS web-build

RUN apt-get update \
    && apt-get install --no-install-recommends -y build-essential ca-certificates curl libssl-dev pkg-config \
    && rm -rf /var/lib/apt/lists/* \
    && curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain 1.88.0

ENV PATH="/root/.cargo/bin:${PATH}"

RUN rustup target add wasm32-unknown-unknown \
    && cargo install wasm-pack --version 0.13.1 --locked

WORKDIR /src
COPY Cargo.toml Cargo.lock LICENSE ./
COPY crates/ crates/
COPY examples/ examples/
COPY web/ web/
RUN sh web/build.sh

FROM python:3.11.18-slim-bookworm@sha256:a83ba18aef0f865b1a4a611c2cf460fe446f2795673bbb9fe43be570439e3b39

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PLAYWRIGHT_BROWSERS_PATH=/ms-playwright

RUN apt-get update \
    && apt-get install --no-install-recommends -y ca-certificates libatomic1 libgcc-s1 libstdc++6 procps r-base-core r-cran-ggplot2 \
    && rm -rf /var/lib/apt/lists/*

RUN python -m pip install --no-cache-dir snakemake==9.14.3

COPY benchmarks/requirements.txt /opt/benchmark/requirements.txt
RUN python -m venv /opt/python-clinker \
    && /opt/python-clinker/bin/pip install --no-cache-dir \
        --requirement /opt/benchmark/requirements.txt

COPY --from=rust-build /opt/clinker-cli /opt/clinker-cli
COPY --from=web-build /usr/local/bin/node /usr/local/bin/node
COPY --from=web-build /src/web/dist /opt/clinker-web
COPY --from=web-build /src/web/node_modules/playwright /opt/benchmark/node_modules/playwright
COPY --from=web-build /src/web/node_modules/playwright-core /opt/benchmark/node_modules/playwright-core
COPY benchmarks/scripts/ /opt/benchmark/
COPY benchmarks/Snakefile benchmarks/fixtures.tsv /opt/benchmark/workflow/
COPY benchmarks/fixtures/ /opt/benchmark/workflow/fixtures/

RUN ln -s /opt/clinker-cli/bin/clinker /usr/local/bin/clinker-rust \
    && node /opt/benchmark/node_modules/playwright/cli.js install --with-deps --only-shell chromium \
    && rm -rf /var/lib/apt/lists/*

RUN { sha256sum /opt/clinker-cli/bin/clinker /opt/clinker-web/build.json \
      /opt/benchmark/requirements.txt /opt/benchmark/*.py /opt/benchmark/*.mjs \
      /opt/benchmark/*.R /opt/benchmark/workflow/Snakefile; \
      /opt/python-clinker/bin/pip freeze; python -m pip freeze; node --version; dpkg-query -W; \
    } | sha256sum | cut -d ' ' -f1 > /opt/benchmark/build-fingerprint

LABEL org.opencontainers.image.title="Clinker benchmark tools" \
      org.opencontainers.image.description="Pinned Clinker CLIs, browser app, and benchmark runners"

WORKDIR /work
ENTRYPOINT ["python", "/opt/benchmark/run.py"]
