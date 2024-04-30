def process_file(input_file, output_file):
    with open(input_file, 'r') as f_input:
        with open(output_file, 'w') as f_output:
            for line in f_input:
                # Split the line into elements
                elements = line.strip().split()
                
                # Enclose each element in double quotes and join with commas
                csv_line = ','.join(['"{}",'.format(elem) for elem in elements])
                
                # Write the processed line to the output file
                f_output.write(csv_line)

def main():
    input_file = input("Enter the path of the input file: ")
    output_file = input("Enter the path of the output file: ")
    
    process_file(input_file, output_file)
    
    print("File processing complete.")

if __name__ == "__main__":
    main()