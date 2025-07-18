#!/usr/bin/env python3
'''Dummy TCP serve to test communication speed.'''
import socket
import threading
import numpy as np

SERVER_IP = "127.0.0.1"
SERVER_PORT = 32409
N_LISTEN = 5

PRE_LEN = 100_000_000
SEED = 0

def send_dummy(client, length, data):
    '''Send dummy data to the client.

    Parameters
    ----------
    client : socket
        Client socket object.
    length : int
        Length of the data to be sent.
    '''
    client.send(data[:length])
    client.close()


class DummyServer:
    '''TCP server that sends dummy data.'''
    def __init__(self, data, ip_address=SERVER_IP, port_num=SERVER_PORT):
        self._sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self._sock.bind((ip_address, port_num))
        self._data = data

    def listen(self, n_listen=N_LISTEN):
        '''Wait connection.'''
        self._sock.listen(n_listen)

        while True:
            client, _ = self._sock.accept()
            req = client.recv(1024)
            try:
                length = int(req)
                assert 0 < length < PRE_LEN
            except AssertionError:
                client.close()
                continue
            except ValueError:
                client.close()
                continue

            client_handler = threading.Thread(target=send_dummy,
                                              args=(client,int(req), self._data))
            client_handler.start()


def main():
    '''Main function'''
    rand_state = np.random.RandomState(SEED)
    data = rand_state.randint(0, 256, size=100_000_000, dtype=np.uint8)
    data = bytes(data)

    server = DummyServer(data)
    server.listen()


if __name__ == '__main__':
    main()
