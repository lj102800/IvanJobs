# Echo server program
import socket

HOST = "127.0.0.1"                 # Symbolic name meaning all available interfaces
PORT = 50007              # Arbitrary non-privileged port

if __name__ == "__main__":
	s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	s.bind((HOST, PORT))
	print "bind (%s, %d)" % (HOST, PORT)
	s.listen(1)
	while True:
		conn, addr = s.accept()
		print "Connected by ", addr
		while True:
			data = conn.recv(1024)
			print "recieved data:%s" % repr(data)
			if not data: break
			conn.send(data)
		conn.close()

