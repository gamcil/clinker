import http.server
import socketserver
import webbrowser
import json
import shutil
import pathlib
import logging

from functools import partial


LOG = logging.getLogger(__name__)


class CustomHandler(http.server.BaseHTTPRequestHandler):
    """Handler for serving cblaster plots."""

    def __init__(self, data, *args, **kwargs):
        self._data = data
        self._dir = pathlib.Path(__file__).resolve().parent.parent / "clinker" / "plot"
        super().__init__(*args, **kwargs)

    def copy_file(self, source):
        shutil.copyfileobj(source, self.wfile)

    def send_headers(self, mime):
        self.send_response(200)
        self.send_header("Content-Type", mime)
        self.end_headers()

    def log_message(self, format, *args):
        """Suppresses logging messages on every request."""
        return

    def do_GET(self):
        """Serves each component of the cblaster plot."""
        if self.path == "/data.json":
            self.send_headers("text/json")
            self.wfile.write(json.dumps(self._data).encode())
            return
        path, mime = None, None
        if self.path == "/":
            path, mime = self._dir / "index.html", "text/html"
        elif self.path == "/style.css":
            path, mime = self._dir / "style.css", "text/css"
        elif self.path == "/d3.min.js":
            path, mime = self._dir / "d3.min.js", "text/javascript"
        elif self.path == "/clustermap.min.js":
            path, mime = self._dir / "clustermap.min.js", "text/javascript"
        elif self.path == "/clinker.js":
            path, mime = self._dir / "clinker.js", "text/javascript"
        if not path:
            return
        with path.open("rb") as fp:
            self.send_headers(mime)
            self.copy_file(fp)


def serve_html(data):
    """Serve a synthaser plot using the socketserver module."""
    handler = partial(CustomHandler, data)

    # Instantiate a new server, bind to any open port
    with socketserver.TCPServer(("localhost", 0), handler) as httpd:

        # Automatically open web browser to bound address
        address, port = httpd.server_address
        url = f"http://{address}:{port}/"
        webbrowser.open(url)

        # Start serving the plot; shutdown on a keyboard interrupt
        try:
            LOG.info(f"Serving clinker plot at {url} (Ctrl+C to stop).")
            httpd.serve_forever()
        except KeyboardInterrupt:
            httpd.shutdown()


def save_html(data, output):
    """Generates a static HTML file with all visualisation code."""

    directory = pathlib.Path(__file__).resolve().parent.parent / "clinker" / "plot"

    with (directory / "index.html").open() as fp:
        html = fp.read()

    css_string = '<link rel="stylesheet" href="style.css"></link>'
    d3_string = '<script src="d3.min.js"></script>'
    cl_string = '<script src="clinker.js"></script>'
    cm_string = '<script src="clustermap.min.js"></script>'

    with (directory / "style.css").open() as fp:
        css = fp.read()
        html = html.replace(css_string, f"<style>{css}</style>")

    with (directory / "d3.min.js").open() as fp:
        d3 = fp.read()
        html = html.replace(d3_string, f"<script>{d3}</script>")

    with (directory / "clustermap.min.js").open() as fp:
        cm = fp.read()
        html = html.replace(cm_string, f"<script>{cm}</script>")

    with (directory / "clinker.js").open() as fp:
        cl = f"const data={json.dumps(data)};" + fp.read()
        html = html.replace(cl_string, f"<script>{cl}</script>")

    with open(output, "w") as fp:
        fp.write(html)


def plot_clusters(clusters, output=None, use_file_order=False):
    """Generates clinker plot from a collection of Synthase objects."""
    data = clusters.to_data(use_file_order=use_file_order)
    plot_data(data, output=output)


def plot_data(data, output=None):
    if output:
        save_html(data, output)
    else:
        serve_html(data)
