import argparse, socket

parser = argparse.ArgumentParser(prog = 'recv', description = 'Send file.')
parser.add_argument('-b', '--buffer', type = int, dest = 'buffer', default = 1024)
parser.add_argument('-p', '--port', type = int, dest = 'port', default = 7878)
parser.add_argument('-v', '--verbose', action = 'store_true')
parser.add_argument('input', type = str)
parser.add_argument('address', type = str)
args = parser.parse_args()

content = open(args.input, 'rb')
address = (args.address, args.port)

s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect(address)

if args.verbose:
	print('Connected to %s:%d' % address)

sent = 0
total = 0
chunk = content.read(args.buffer)

while True:
	sent = s.send(chunk)
	chunk = chunk[sent:]
	total += sent

	part = content.read(sent)
	
	if part:
		chunk += part

	if args.verbose:
		print('\rSent %d bytes' % total, end = '')

	if not chunk:
		break
	
if args.verbose:
	print('\nDone')
