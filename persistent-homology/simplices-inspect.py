# Python script to find duplicate lines in a file

def find_duplicate_lines(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    line_count = {}
    for line in lines:
        line = line.strip()
        if line in line_count:
            line_count[line] += 1
        else:
            line_count[line] = 1

    duplicates = {line: count for line, count in line_count.items() if count > 1}

    if duplicates:
        for line, count in duplicates.items():
            print(f"{line} (appears {count} times)")
        print(f"{len(duplicates)} duplicate lines found:")
    else:
        print("No duplicate lines found.")

if __name__ == "__main__":
    find_duplicate_lines('st-visualizer/tmp.txt')