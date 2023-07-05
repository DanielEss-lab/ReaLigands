class HeaderReplacer:
    def __init__(self, filename):
        self.filename = filename

    def replace_header(self, indices):
        # instead of searching for the header line from the original tmqm files (as done below), this code
        # just replaces the second line of the file with the indices of the connecting atoms
        with open(self.filename, 'r') as f:
           file_data = f.readlines()
           file_data[1] = ','.join([str(i) for i in indices]) + '\n'

        with open(self.filename, 'w') as file:
            file.writelines(file_data)

        # file_data = file.read()
        # header = re.findall("CSD_code(.*)", file_data)[0]

        # replace_string = "CSD_code" + header
        # # new_string = str(index)
        # new_string = ','.join([str(i) for i in indices])

        # file_data = file_data.replace(replace_string, new_string)

        
