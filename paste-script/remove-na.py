import csv

def filter_rows(input_file, output_file):
    with open(input_file, 'r', newline='') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        header = next(reader)
        writer.writerow(header)

        for row in reader:
            last_columns = row[-10:]
            if not all(value == "NA" for value in last_columns):
                row = list(map(lambda x: x.replace("NA", "0"), row))
                writer.writerow(row)

if __name__ == "__main__":
    input_file_path = "old.tsv"
    output_file_path = "new.tsv"
    filter_rows(input_file_path, output_file_path)

