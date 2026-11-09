FROM rust:1.88.0-slim-bookworm@sha256:38bc5a86d998772d4aec2348656ed21438d20fcdce2795b56ca434cf21430d89 AS rust-build

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
