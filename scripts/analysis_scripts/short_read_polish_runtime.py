import sys
import re
import csv
def parse_directory_error(error_message):
    # Split the error message by single quotes
    split_message = error_message.split("‘")

    # Extract the directory name from the split message
    directory_name = split_message[1].split("’")[0].split("/")[-1]

    return directory_name


if len(sys.argv) != 2:
    print("Usage: python program.py input_file")
    sys.exit(1)

filename = sys.argv[1]

with open(filename, "r") as file:
    text = file.read()

first_row = text.strip().split("\n")[0]
method_name = first_row.split(":")[1].strip()
if "cannot create directory" in method_name:
    method_name = parse_directory_error(method_name)
    #print(parsed_directory_name)



print(method_name)


pattern = r'real\s+(\d+m\d+\.\d+s)\s*[\s\S]*?Done - ([\w\s]+) round:\s+(\d+)'

matches = re.findall(pattern, text)
round_times = {}
for match in matches:
    round_type = match[1].strip().lower()
    round_num = int(match[2])
    round_time = match[0]
    round_id = f"{round_type}_{round_num}"
    round_times[round_id] = round_time

output_file = method_name+"_short_read_runtime.csv"

with open(output_file, "w", newline="") as file:
    writer = csv.writer(file)
    writer.writerow(["Method", "Round", "Time"])
    for round_id, round_time in round_times.items():
        writer.writerow([method_name,round_id, round_time])

print(f"Output written to {output_file}")
