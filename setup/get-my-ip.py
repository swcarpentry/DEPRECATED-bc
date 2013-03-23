#!/usr/bin/env python

"""Get your public IP address from a UDP socket connection
"""

import socket as _socket


def get_my_ip(host, port=80):
    s = _socket.socket(_socket.AF_INET, _socket.SOCK_DGRAM)
    try:
        s.connect((host, port))
        return s.getsockname()[0]
    finally:
        s.close()


if __name__ == '__main__':
    import argparse as _argparse

    parser = _argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        'host', default='software-carpentry.org', nargs='?')
    parser.add_argument(
        'port', default=80, type=int, nargs='?')

    args = parser.parse_args()

    print(get_my_ip(host=args.host, port=args.port))
