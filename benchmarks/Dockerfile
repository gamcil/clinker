FROM rust:1.87.0-slim-bookworm@sha256:437507c3e719e4f968033b88d851ffa9f5aceeb2dcc2482cc6cb7647811a55eb AS rust-build

ARG CLINKER_CLI_VERSION=0.1.0
RUN cargo install clinker-cli \
    --version "${CLINKER_CLI_VERSION}" \
    --locked \
    --root /opt/clinker-cli

FROM python:3.11.18-slim-bookworm@sha256:a83ba18aef0f865b1a4a611c2cf460fe446f2795673bbb9fe43be570439e3b39

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

COPY requirements.txt /opt/benchmark/requirements.txt
RUN python -m venv /opt/python-clinker \
    && /opt/python-clinker/bin/pip install --no-cache-dir \
        --requirement /opt/benchmark/requirements.txt

COPY --from=rust-build /opt/clinker-cli /opt/clinker-cli
COPY scripts/ /opt/benchmark/

RUN ln -s /opt/clinker-cli/bin/clinker /usr/local/bin/clinker-rust

LABEL org.opencontainers.image.title="Clinker benchmark tools" \
      org.opencontainers.image.description="Pinned Python and Rust Clinker CLIs and benchmark measurement scripts"
