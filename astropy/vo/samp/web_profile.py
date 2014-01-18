from ...extern.six.moves.urllib.parse import parse_qs
from ...extern.six.moves.urllib.request import urlopen
from ...extern.six.moves import input

from .standard_profile import SAMPSimpleXMLRPCRequestHandler, ThreadingXMLRPCServer

__all__ = []


class WebProfileRequestHandler(SAMPSimpleXMLRPCRequestHandler):
    """
    Handler of XMLRPC requests performed through the Web Profile.
    """

    def _send_CORS_header(self):

        if not self.headers.get('Origin') is None:

            method = self.headers.get('Access-Control-Request-Method')
            if method and self.command == "OPTIONS":
                # Preflight method
                self.send_header('Content-Length', '0')
                self.send_header('Access-Control-Allow-Origin', self.headers.get('Origin'))
                self.send_header('Access-Control-Allow-Methods', method)
                self.send_header('Access-Control-Allow-Headers', 'Content-Type')
                self.send_header('Access-Control-Allow-Credentials', 'true')
            else:
                # Simple method
                self.send_header('Access-Control-Allow-Origin', self.headers.get('Origin'))
                self.send_header('Access-Control-Allow-Headers', 'Content-Type')
                self.send_header('Access-Control-Allow-Credentials', 'true')

    def end_headers(self):
        self._send_CORS_header()
        SAMPSimpleXMLRPCRequestHandler.end_headers(self)

    def _serve_cross_domain_xml(self):

        cross_domain = False

        if self.path == "/crossdomain.xml":
            # Adobe standard
            response = """<?xml version='1.0'?>
<!DOCTYPE cross-domain-policy SYSTEM "http://www.adobe.com/xml/dtds/cross-domain-policy.dtd">
<cross-domain-policy>
  <site-control permitted-cross-domain-policies="all"/>
  <allow-access-from domain="*"/>
  <allow-http-request-headers-from domain="*" headers="*"/>
</cross-domain-policy>"""

            self.send_response(200, 'OK')
            self.send_header('Content-Type', 'text/x-cross-domain-policy')
            self.send_header("Content-Length", str(len(response)))
            self.end_headers()
            self.wfile.write(response)
            self.wfile.flush()
            cross_domain = True

        elif self.path == "/clientaccesspolicy.xml":
            # Microsoft standard
            response = """<?xml version='1.0'?>
<access-policy>
  <cross-domain-access>
    <policy>
      <allow-from>
        <domain uri="*"/>
      </allow-from>
      <grant-to>
        <resource path="/" include-subpaths="true"/>
      </grant-to>
    </policy>
  </cross-domain-access>
</access-policy>"""

            self.send_response(200, 'OK')
            self.send_header('Content-Type', 'text/xml')
            self.send_header("Content-Length", str(len(response)))
            self.end_headers()
            self.wfile.write(response)
            self.wfile.flush()
            cross_domain = True

        return cross_domain

    def do_POST(self):
        if self._serve_cross_domain_xml():
            return

        return SAMPSimpleXMLRPCRequestHandler.do_POST(self)

    def do_HEAD(self):

        if not self.is_http_path_valid():
            self.report_404()
            return

        if self._serve_cross_domain_xml():
            return

    def do_OPTIONS(self):

        self.send_response(200, 'OK')
        self.end_headers()

    def do_GET(self):

        if not self.is_http_path_valid():
            self.report_404()
            return

        split_path = self.path.split('?')

        if split_path[0] in ['/translator/%s' % clid for clid in self.server.clients]:
            # Request of a file proxying
            urlpath = parse_qs(split_path[1])
            try:
                proxyfile = urlopen(urlpath["ref"][0])
                self.send_response(200, 'OK')
                self.end_headers()
                self.wfile.write(proxyfile.read())
                proxyfile.close()
            except:
                self.report_404()
                return

        if self._serve_cross_domain_xml():
            return

    def is_http_path_valid(self):

        valid_paths = ["/clientaccesspolicy.xml", "/crossdomain.xml"] + ['/translator/%s' % clid for clid in self.server.clients]
        return self.path.split('?')[0] in valid_paths


class WebProfileXMLRPCServer(ThreadingXMLRPCServer):
    """
    XMLRPC server supporting the SAMP Web Profile.
    """

    def __init__(self, addr, log=None, requestHandler=WebProfileRequestHandler,
                 logRequests=True, allow_none=True, encoding=None):

        self.clients = []
        ThreadingXMLRPCServer.__init__(self, addr, log, requestHandler,
                                       logRequests, allow_none, encoding)

    def add_client(self, client_id):
        self.clients.append(client_id)

    def remove_client(self, client_id):
        try:
            self.clients.remove(client_id)
        except:
            pass


def web_profile_text_dialog(request, queue):

    samp_name = "unknown"

    if isinstance(request[0], str):
        # To support the old protocol version
        samp_name = request[0]
    else:
        samp_name = request[0]["samp.name"]

    text = \
        """A Web application which declares to be

Name: %s
Origin: %s

is requesting to be registered with the SAMP Hub.
Pay attention that if you permit its registration, such
application will acquire all current user privileges, like
file read/write.

Do you give your consent? [yes|no]""" % (samp_name, request[2])

    print(text)
    answer = input(">>> ")
    queue.put(answer.lower() in ["yes", "y"])
