from preprocessing import *
import json
import sys

def main(arg):
	params = json.loads(arg)
	task = PreprocessChunkTask(**params)
	task.execute()

if __name__ == "__main__":
	main(sys.argv[1])