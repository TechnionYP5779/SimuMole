import csv


def print_clean_status(path):
    print("  0: " + read_clean_status(path, 'SimulationForm0_LoadPdb'))
    print("  1: " + read_clean_status(path, 'SimulationForm1_DetermineRelativePosition'))
    print("  2: " + read_clean_status(path, 'SimulationForm2_SimulationParameters'))


def init_clean_status(path, cleaned_data={}):
    with open(path, mode='w', newline='') as csv_file:
        # status: {pass, fail, unknown}
        fieldnames = ['form_name', 'status']
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames, delimiter=',')
        writer.writeheader()

        # ................. form 0 ................. #
        writer.writerow({'form_name': 'SimulationForm0_LoadPdb', 'status': 'unknown'})

        writer.writerow({'form_name': 'num_of_proteins', 'status': str(cleaned_data.get('num_of_proteins', '?'))})
        writer.writerow({'form_name': 'first_pdb_type', 'status': str(cleaned_data.get('first_pdb_type', '?'))})
        writer.writerow({'form_name': 'first_pdb_id', 'status': str(cleaned_data.get('first_pdb_id', '?'))})
        writer.writerow({'form_name': 'first_pdb_file', 'status': str(cleaned_data.get('first_pdb_file', '?'))})
        writer.writerow({'form_name': 'second_pdb_type', 'status': str(cleaned_data.get('second_pdb_type', '?'))})
        writer.writerow({'form_name': 'second_pdb_id', 'status': str(cleaned_data.get('second_pdb_id', '?'))})
        writer.writerow({'form_name': 'second_pdb_file', 'status': str(cleaned_data.get('second_pdb_file', '?'))})

        writer.writerow({'form_name': 'clean_first_pdb_type', 'status': 'unknown'})
        writer.writerow({'form_name': 'clean_first_pdb_id', 'status': 'unknown'})
        writer.writerow({'form_name': 'clean_first_pdb_file', 'status': 'unknown'})
        writer.writerow({'form_name': 'clean_second_pdb_type', 'status': 'unknown'})
        writer.writerow({'form_name': 'clean_second_pdb_id', 'status': 'unknown'})
        writer.writerow({'form_name': 'clean_second_pdb_file', 'status': 'unknown'})

        # ................. form 1 ................. #
        writer.writerow({'form_name': 'SimulationForm1_DetermineRelativePosition', 'status': 'unknown'})

        # ................. form 2 ................. #
        writer.writerow({'form_name': 'SimulationForm2_SimulationParameters', 'status': 'unknown'})


def read_clean_status(path, required_form_name):
    with open(path) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                pass
            else:
                form_name = row[0]
                status = row[1]
                if required_form_name == form_name:
                    return status
            line_count += 1


def write_clean_status(path, required_form_name, required_status):
    with open(path) as in_csv_file:
        reader = csv.reader(in_csv_file.readlines())

    with open(path, 'w', newline='') as out_in_csv_file:
        writer = csv.writer(out_in_csv_file)
        for line in reader:
            if line[0] == required_form_name:
                writer.writerow([required_form_name, required_status])
                break
            else:
                writer.writerow(line)
        writer.writerows(reader)


def same_form__SimulationForm0_LoadPdb(path, cleaned_data):
    return (read_clean_status(path, 'num_of_proteins') == str(cleaned_data.get('num_of_proteins', '?'))) and \
           (read_clean_status(path, 'first_pdb_type') == str(cleaned_data.get('first_pdb_type', '?'))) and \
           (read_clean_status(path, 'first_pdb_id') == str(cleaned_data.get('first_pdb_id', '?'))) and \
           (read_clean_status(path, 'first_pdb_file') == str(cleaned_data.get('first_pdb_file', '?'))) and \
           (read_clean_status(path, 'second_pdb_type') == str(cleaned_data.get('second_pdb_type', '?'))) and \
           (read_clean_status(path, 'second_pdb_id') == str(cleaned_data.get('second_pdb_id', '?'))) and \
           (read_clean_status(path, 'second_pdb_file') == str(cleaned_data.get('second_pdb_file', '?')))


def complete_cleaning__SimulationForm0_LoadPdb(path, cleaned_data):
    num_of_proteins = cleaned_data['num_of_proteins']
    first_pdb_type = cleaned_data['first_pdb_type']
    second_pdb_type = cleaned_data['second_pdb_type']

    clean_first_pdb_type = read_clean_status(path, "clean_first_pdb_type") == "pass"
    clean_first_pdb_id = read_clean_status(path, "clean_first_pdb_id") == "pass"
    clean_first_pdb_file = read_clean_status(path, "clean_first_pdb_file") == "pass"
    clean_second_pdb_type = read_clean_status(path, "clean_second_pdb_type") == "pass"
    clean_second_pdb_id = read_clean_status(path, "clean_second_pdb_id") == "pass"
    clean_second_pdb_file = read_clean_status(path, "clean_second_pdb_file") == "pass"

    first_pdb_status = (first_pdb_type == 'by_id' and clean_first_pdb_id) or \
                       (first_pdb_type == 'by_file' and clean_first_pdb_file)
    second_pdb_status = (second_pdb_type == 'by_id' and clean_second_pdb_id) or \
                        (second_pdb_type == 'by_file' and clean_second_pdb_file)

    return clean_first_pdb_type and clean_second_pdb_type and \
           ((num_of_proteins == '1' and first_pdb_status)
            or
            (num_of_proteins == '2' and first_pdb_status and second_pdb_status))
