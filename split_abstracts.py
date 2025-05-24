import os
import re
import shutil

def clear_directory(directory):
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.remove(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print(f"Failed to delete {file_path}. Reason: {e}")

def split_file(input_file, delimiter_pattern, journal_name):
    base_directory = "abstracts"
    output_directory = f"{base_directory}/{input_file}_output"
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    with open(input_file, 'r') as file:
        content = file.read()

    # Adjust the regular expression to not capture the delimiter
    sections = re.split(delimiter_pattern, content)

    # Initialize a regex to identify minimal valid content
    minimal_content_regex = re.compile(r'\w+')

    for i, section in enumerate(sections, 1):
        section = section.strip()
        if section and minimal_content_regex.search(section):  # Check if the section contains word characters
            output_file = f"{output_directory}/output_section_{i - 1}.txt"
            content_to_write = section
            if i != 1:
                content_to_write = f"Journal={journal_name}, " + section  # Prepend journal name to each section except the first one
            # Further checking against only containing the journal name and delimiter
            if minimal_content_regex.search(content_to_write[len(f"Journal={journal_name}, "):]):
                with open(output_file, 'w') as file:
                    file.write(content_to_write)
                print(f"Section {i} from {input_file} written to {output_file}")

# Clear the abstracts directory at the start
abstracts_directory = "abstracts"
if os.path.exists(abstracts_directory):
    clear_directory(abstracts_directory)
else:
    os.makedirs(abstracts_directory)

# Example usage:

input_files = [
    ('abstract-AmericanJo-set.txt', r'\d{1,5}\. Am J Epidemiol\.', 'Am J Epidemiol'),
    ('abstract-Epidemiolo-set.txt', r'\d{1,5}\. Epidemiology\.', 'Epidemiology'),
    ('abstract-EuropeanJo-set.txt', r'\d{1,5}\. Eur J Epidemiol\.', 'Eur J Epidemiol'),
    ('abstract-Internatio-set.txt', r'\d{1,5}\. Int J Epidemiol\.', 'Int J Epidemiol')
]

for input_file, delimiter_pattern, journal_name in input_files:
    split_file(input_file, delimiter_pattern, journal_name)
