import subprocess

with open("foobar.txt", "w") as outh:
	command = ["grep", "-wF", "MSTRG.2.1", "/scratch/user/uqktiwar/NCCSOP/merged/soxPN.f0.01.gtf"]

	try:
		output = subprocess.check_output(command, stderr=subprocess.STDOUT, text=True)
		print(output)
		outh.write(output)
	except subprocess.CalledProcessError as e:
		print("Error:", e.output)

	command = ["grep", "-wF", "MSTRG.2.3", "/scratch/user/uqktiwar/NCCSOP/merged/soxPN.f0.01.gtf"]

	try:
		output = subprocess.check_output(command, stderr=subprocess.STDOUT, text=True)
		print(output)
		outh.write(output)
	except subprocess.CalledProcessError as e:
		print("Error:", e.output)



