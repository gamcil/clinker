use std::{io, io::Read, process::Command};

use clinker_core::{AnalysisSession, ClinkerProject, ClusterPairSimilarity};
use tiny_http::{Header, Method, Response, Server, StatusCode};

include!(concat!(env!("OUT_DIR"), "/viewer-assets.rs"));

/// Serve only the embedded viewer and explicit JSON endpoints, never paths on
/// disk. A fresh loopback port and capability path identify this running figure.
pub fn serve(
    project: &ClinkerProject,
    similarities: &[Vec<ClusterPairSimilarity>],
    no_browser: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let server =
        Server::http("127.0.0.1:0").map_err(|error| io::Error::other(error.to_string()))?;
    let address = server.server_addr().to_ip().expect("TCP listener");
    let mut random = [0_u8; 16];
    getrandom::getrandom(&mut random).map_err(|error| io::Error::other(error.to_string()))?;
    let prefix = format!(
        "/{}/",
        random
            .iter()
            .map(|byte| format!("{byte:02x}"))
            .collect::<String>()
    );
    let origin = format!("http://{address}");
    let url = format!("{origin}{prefix}");
    let project_json = project.to_json()?;
    let similarity = serde_json::to_vec(similarities)?;
    println!("Viewer: {url}\nPress Ctrl+C to stop the local viewer.");
    if !no_browser && let Err(error) = open_browser(&url) {
        eprintln!("Could not open a browser ({error}). Open the URL above manually.");
    }

    for mut request in server.incoming_requests() {
        let host_matches = request.headers().iter().any(|header| {
            header.field.equiv("Host") && header.value.as_str() == address.to_string()
        });
        let origin_matches = request
            .headers()
            .iter()
            .all(|header| !header.field.equiv("Origin") || header.value.as_str() == origin);
        let path = request.url().split('?').next().unwrap_or("");
        let Some(path) = path
            .strip_prefix(&prefix)
            .filter(|_| host_matches && origin_matches)
        else {
            let _ = request.respond(Response::from_string("Not found").with_status_code(404));
            continue;
        };
        // Own the route before borrowing the request body below.
        let path = path.to_owned();
        let result = match (request.method(), path.as_str()) {
            (Method::Get | Method::Head, "api/project") => {
                Ok((project_json.clone(), "application/json"))
            }
            (Method::Get | Method::Head, "api/similarity") => {
                Ok((similarity.clone(), "application/json"))
            }
            (Method::Post, "api/similarity") => read_similarity(&mut request),
            (Method::Get | Method::Head, path) => {
                let name = if path.is_empty() { "index.html" } else { path };
                ASSETS
                    .iter()
                    .find(|(key, _)| *key == name)
                    .map(|(_, bytes)| (bytes.to_vec(), content_type(name)))
                    .ok_or((404, "Not found".to_owned()))
            }
            _ => Err((405, "Method not allowed".to_owned())),
        };
        let (code, bytes, mime) = match result {
            Ok((bytes, mime)) => (200, bytes, mime),
            Err((code, message)) => (code, message.into_bytes(), "text/plain; charset=utf-8"),
        };
        let response = Response::from_data(bytes)
            .with_status_code(StatusCode(code))
            .with_header(Header::from_bytes("Content-Type", mime).unwrap())
            .with_header(Header::from_bytes("Cache-Control", "no-store").unwrap());
        if let Err(error) = request.respond(response) {
            eprintln!("Viewer connection closed: {error}");
        }
    }
    Ok(())
}

fn read_similarity(
    request: &mut tiny_http::Request,
) -> Result<(Vec<u8>, &'static str), (u16, String)> {
    const MAX_BYTES: u64 = 256 * 1024 * 1024;
    let mut bytes = Vec::new();
    request
        .as_reader()
        .take(MAX_BYTES + 1)
        .read_to_end(&mut bytes)
        .map_err(|error| (400, error.to_string()))?;
    if bytes.len() as u64 > MAX_BYTES {
        return Err((
            413,
            "Analysis session exceeds the 256 MiB request limit".to_owned(),
        ));
    }
    let session = AnalysisSession::from_json(&bytes).map_err(|error| (400, error.to_string()))?;
    let matrix = serde_json::to_vec(&session.analysis.cluster_similarity_matrix())
        .map_err(|error| (400, error.to_string()))?;
    Ok((matrix, "application/json"))
}

fn content_type(name: &str) -> &'static str {
    match name.rsplit('.').next().unwrap_or("") {
        "html" => "text/html; charset=utf-8",
        "js" | "mjs" => "text/javascript; charset=utf-8",
        "css" => "text/css; charset=utf-8",
        "json" => "application/json",
        _ => "text/plain; charset=utf-8",
    }
}

fn open_browser(url: &str) -> io::Result<()> {
    #[cfg(target_os = "macos")]
    let status = Command::new("open").arg(url).status()?;
    #[cfg(target_os = "windows")]
    let status = Command::new("cmd")
        .args(["/C", "start", "", url])
        .status()?;
    #[cfg(not(any(target_os = "macos", target_os = "windows")))]
    let status = Command::new("xdg-open").arg(url).status()?;
    if status.success() {
        Ok(())
    } else {
        Err(io::Error::other("browser launcher failed"))
    }
}
