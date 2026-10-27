//! Self-contained renderer-only HTML export. The full editor uses the local viewer.

pub fn html(plot_json: &str) -> String {
    const D3: &str = include_str!("../assets/viewer/vendor/d3/d3.min.js");
    const CLUSTERMAP: &str = include_str!("../assets/viewer/vendor/clustermap/clustermap.min.js");
    let safe_json = plot_json
        .replace('&', "\\u0026")
        .replace('<', "\\u003c")
        .replace('>', "\\u003e")
        .replace('\u{2028}', "\\u2028")
        .replace('\u{2029}', "\\u2029");
    let mut html = String::with_capacity(D3.len() + CLUSTERMAP.len() + safe_json.len() + 1_024);
    html.push_str(
        "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">\
         <meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">\
         <title>clinker plot</title><style>html,body{width:100%;height:100%;margin:0;overflow:hidden;font-family:system-ui,sans-serif}\
         #plot{position:relative;width:100%;height:100%;overflow:hidden}svg.clusterMap{display:block}</style></head><body>\
         <div id=\"plot\"></div><script>",
    );
    html.push_str(D3);
    html.push_str("</script><script>");
    html.push_str(CLUSTERMAP);
    html.push_str("</script><script>const data=");
    html.push_str(&safe_json);
    html.push_str(
        ";const chart=ClusterMap.ClusterMap();d3.select('#plot').datum(data).call(chart);</script></body></html>",
    );
    html
}
