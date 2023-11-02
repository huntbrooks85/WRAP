"""
Simple to use IPv6 compatible Internet connection functions
"""

from __future__ import print_function

__author__ = 'Phil Budne <phil@ultimate.com>'
__version__ = '0.2'
__revision__ = '$Id: connect.py,v 1.3 2020/06/14 04:32:43 phil Exp $'

import socket

def _conn(host, port, req_fam, req_type):
    """worker function"""
    error = None
    # getaddrinfo will raise exception if no matches
    for (fam, type, proto, cname, addr) in socket.getaddrinfo(host, port, req_fam, req_type):
        try:
            s = socket.socket(fam, type, proto)
            s.connect(addr)
            return s
        except socket.error as e:
            error = e
    # pass up last error
    raise error

def tcp(host, port):
    """return connected TCP socket using either IPv4 or IPv6"""
    return _conn(host, port, socket.AF_UNSPEC, socket.SOCK_STREAM)

def udp(host, port):
    """return connected UDP socket using either IPv4 or IPv6"""
    return _conn(host, port, socket.AF_UNSPEC, socket.SOCK_DGRAM)

def tcp4(host, port):
    """return connected TCP socket using IPv4"""
    return _conn(host, port, socket.AF_INET, socket.SOCK_STREAM)

def udp4(host, port):
    """return connected UDP socket using IPv4"""
    return _conn(host, port, socket.AF_INET, socket.SOCK_DGRAM)

def tcp6(host, port):
    """return connected TCP socket using IPv6"""
    return _conn(host, port, socket.AF_INET6, socket.SOCK_STREAM)

def udp6(host, port):
    """return connected UDP socket using IPv6"""
    return _conn(host, port, socket.AF_INET6, socket.SOCK_DGRAM)

def getfqdn4(name=None):
    """return (a) IPv4 FQDN (Fully Qualified Domain Name)
       if name is not given, returns local hostname"""
    if name is None:
        return socket.getfqdn()
    return socket.getfqdn(name)

def _getfqdn(rfam, n=None):
    """INTERNAL: return FQDN given address family (and name)
       may raise socket.gaierror"""
    p = 80
    if n is None:
        n = socket.gethostname()
    rtype = socket.SOCK_STREAM
    for (fam, type, proto, cname, addr) in socket.getaddrinfo(n,p,rfam,rtype):
        nn = socket.getnameinfo(addr, 0)[0]
        if nn != addr[0]:
            return nn
    return None

def getfqdn6(name=None):
    """return (a) IPv6 FQDN (Fully Qualified Domain Name), given name
       if name is not given, returns local hostname
       may raise socket.gaierror"""
    return _getfqdn(socket.AF_INET6, name)

def getfqdn(name=None):
    """return (a) local IPv4 or v6 FQDN (Fully Qualified Domain Name)
       if name is not given, returns local hostname
       may raise socket.gaierror"""
    return _getfqdn(socket.AF_UNSPEC, name)

if __name__ == '__main__':

    c = tcp4('localhost', 'www')
    if c:
        print(c.getpeername())

    print("getfqdn4", getfqdn4())
    print("getfqdn", getfqdn())
    print("getfqdn6", getfqdn6())
