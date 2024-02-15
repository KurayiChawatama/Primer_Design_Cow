import pyfiglet
import cowsay
import pyttsx3
import time
import requests
import re
from tabulate import tabulate
from Levenshtein import distance
from io import StringIO

engine = pyttsx3.init()


def main():
    global engine
    welcome_sequence()
    accession_number = input("Enter a Gene Accession Number: ")
    gene_summary_data, gene_id = fetch_gene_summary(accession_number)

    if gene_summary_data and gene_id:
        cow_summary(gene_summary_data, gene_id)
        snp, length = cow_pick_snp(gene_summary_data, gene_id)
        sequence_data = fetch_sequence_data(gene_id)
        # Check if sequence data was fetched
        if sequence_data:
            fetched_successfully = "Gene sequence Data Fetched Successfully. Finding suitable primers"
            display_and_speak(fetched_successfully)
            up_flank, down_flank = collect_flanking_nucleotides(sequence_data, snp, length)
            # Check if the flanks are suitable for primer design
            if up_flank and down_flank:
                forward_primers = []
                reverse_primers = []
                forward_primers, forward_positions = find_primers(up_flank, forward_primers, "forward")
                reverse_primers, reverse_positions = find_primers(down_flank, reverse_primers, "reverse")
                # Adjust positions
                forward_positions = [pos + max(0, snp - 100) for pos in forward_positions]
                reverse_positions = [pos + snp for pos in reverse_positions]
                # Calculate CG percentages
                forward_cg_percentages = []
                reverse_cg_percentages = []
                forward_cg_percentages = calculate_cg_percentage(forward_primers, forward_cg_percentages)
                reverse_cg_percentages = calculate_cg_percentage(reverse_primers, reverse_cg_percentages)
                # Count nucleotides
                forward_nucleotide_count = count_nucleotides(forward_primers)
                reverse_nucleotide_count = count_nucleotides(reverse_primers)
                # Calculate melting temperatures
                forward_meting_temps = []
                reverse_melting_temps = []
                reverse_melting_temps = calculate_melting_temps(reverse_primers, reverse_melting_temps)
                forward_meting_temps = calculate_melting_temps(forward_primers, forward_meting_temps)
                # Print primer tables
                forward_table = print_primer_table(forward_primers, forward_cg_percentages, forward_nucleotide_count,
                                                   forward_meting_temps, "Forward", forward_positions)
                reverse_table = print_primer_table(reverse_primers, reverse_cg_percentages, reverse_nucleotide_count,
                                                   reverse_melting_temps, "Reverse", reverse_positions)
                forward_tables_message = "The forward primer table is displayed below."
                reverse_tables_message = "The reverse primers table is displayed below"
                display_and_speak(forward_tables_message)
                print(forward_table)
                display_and_speak(reverse_tables_message)
                print(reverse_table)

                select_primer_pair(forward_primers, reverse_primers, forward_meting_temps, reverse_melting_temps)
                # Save files
                while True:
                    save_prompt = "Would you like to save the files?"
                    display_and_speak(save_prompt)
                    save_decision = input("Would you like to save the files? (y/n): ")
                    if save_decision.lower() == "y":
                        files_saved = "Files saved."
                        save_all_files(accession_number, snp, forward_table, reverse_table, sequence_data)
                        display_and_speak(files_saved)
                        break
                    elif save_decision.lower() == "n":
                        files_not_saved = "Files not saved."
                        display_and_speak(files_not_saved)
                        break
                    else:
                        invalid_input = "Invalid input. select y or n."
                        display_and_speak(invalid_input)


def display_and_speak(*messages):
    global engine
    message = ''.join(messages)
    cowsay.cow(message)
    engine.say(message)
    engine.runAndWait()


def print_figlet_text(text):
    print(pyfiglet.figlet_format(text))
    time.sleep(1)


def welcome_sequence():
    print_figlet_text("The")
    print_figlet_text("Primer Design")
    print_figlet_text("Cow!")
    welcome = "Welcome to the primer design cow! Please enter a gene accession number to get started."
    display_and_speak(welcome)


def fetch_gene_summary(accession):
    fetching_data_message = "Fetching Gene Sequence Data from NCBI Database"
    display_and_speak(fetching_data_message)
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id={accession}&retmode=json"
    try:
        response = requests.get(url)
        response = response.json()
        try:
            gene_id = response['result']['uids'][0]
            return response, gene_id
        except IndexError:
            failed = f"Failed to fetch gene summary for accession number {accession}. Recheck the accession number."
            display_and_speak(failed)
            raise ValueError()
    except requests.RequestException:
        bad_internet = "Failed to fetch gene summary. Please check your internet connection."
        display_and_speak(bad_internet)
        raise ValueError()


def cow_summary(summary_json, gene_id):
    global engine
    gene_name = summary_json['result'][gene_id]['title']
    date_updated = summary_json['result'][gene_id]['updatedate']
    gen_length = summary_json['result'][gene_id]['slen']
    summary_to_say = f"""The gene is {gene_name}.
    It was last updated on {date_updated}. 
    It is {gen_length} nucleotides long."""
    cowsay.cow(summary_to_say)
    engine.say(summary_to_say)
    engine.runAndWait()


def cow_pick_snp(summary_json, gene_id):
    picker = "Please pick a nucleotide position of interest."
    invalid = "Invalid position. Please enter a position 75 BP within the sequence length for a visible band.(>50BP)"
    display_and_speak(picker)

    gen_length = summary_json['result'][gene_id]['slen']
    position = -1
    while position < 75 or position >= (gen_length - 75):
        try:
            position = int(input("Enter a nucleotide position: "))
        except ValueError:
            display_and_speak(invalid)
            continue
        if position < 75 or position >= (gen_length - 75):
            display_and_speak(invalid)

    position_picked_message = f"You picked the nucleotide position {position}"
    display_and_speak(position_picked_message)
    return position, gen_length


def fetch_sequence_data(gene_id):
    bad_internet = "Failed to fetch gene sequence. Please check your internet connection."
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={gene_id}&rettype=fasta&retmode=text"
    try:
        response = requests.get(url)
        return response.text
    except requests.RequestException:
        display_and_speak(bad_internet)
        raise ValueError()


def collect_flanking_nucleotides(sequence, position, length):
    sequence = sequence.split('\n', 1)[1]
    sequence = sequence.replace('\n', '')

    if position - 100 < 0:
        up_flank = sequence[:position]
    else:
        up_flank = sequence[position - 100:position]

    if position + 100 > length:
        down_flank = sequence[position:]
    else:
        down_flank = sequence[position:position + 100]

    return up_flank, down_flank


def count_nucleotides(primer_list):
    nucleotide_counts = []
    for primer in primer_list:
        nucleotide_count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for nucleotide in primer:
            if nucleotide in nucleotide_count:
                nucleotide_count[nucleotide] += 1
        nucleotide_counts.append(nucleotide_count)
    return nucleotide_counts


def find_primers(flank, primer_list, direction):
    primer = ""
    primer_positions = []
    if direction == "reverse":
        flank = reverse_complement(flank)
    for i in range(len(flank)):
        for length in range(18, 26):  # Loop over possible primer lengths
            if i + length > len(flank):  # Skip if the end of the sequence is reached
                break
            if re.match(r'^[GC][ATCG]{21,23}[GC]$', flank[i:i + length]):
                primer = flank[i:i + length]
                primer_list.append(primer)
                primer_positions.append(i)
    return primer_list, primer_positions


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return "".join(complement[base] for base in reversed(seq))


def calculate_cg_percentage(primer_list, cg_percentages):
    for primer in primer_list:
        cg_count = 0
        for nucleotide in primer:
            if nucleotide == 'C' or nucleotide == 'G':
                cg_count += 1
        cg_percentage = (cg_count / len(primer)) * 100
        cg_percentages.append(cg_percentage)
    return cg_percentages


def calculate_melting_temps(primer_list, melting_temps):
    for primer in primer_list:
        gc_count = 0
        at_count = 0
        for nucleotide in primer:
            if nucleotide == 'G' or nucleotide == 'C':
                gc_count += 1
            if nucleotide == 'A' or nucleotide == 'T':
                at_count += 1
        melting_temp = (gc_count * 4) + (at_count * 2)
        melting_temps.append(melting_temp)
    return melting_temps


def print_primer_table(primer_list, cg_percentages, nucleotide_count, melting_temps, direction, positions):
    table = []
    for i in range(len(primer_list)):
        primer = primer_list[i]
        primer_length = len(primer)  # Calculate the length of the primer
        start_position = positions[i]
        end_position = start_position + primer_length - 1  # Subtract 1 because positions are 0-indexed
        position_range = f"{start_position}-{end_position}"
        if direction == "Reverse":
            primer = reverse_complement(primer)
        table.append([i+1, primer, primer_length, cg_percentages[i], nucleotide_count[i], melting_temps[i],
                      position_range])
    output = StringIO()
    output.write(f"{direction} Primer Table\n")
    output.write(tabulate(table, headers=["Primer Number", "Primer Sequence", "Primer Length", "CG Percentage",
                                          "Nucleotide Count", "Melting Temperature", "Position Range in Gene Sequence"],
                          tablefmt="fancy_grid"))
    output.write("\n")
    return output.getvalue()


def select_primer_pair(forward_primers, reverse_primers, forward_temps, reverse_temps):
    please_select_primer = "Please select a forward and reverse primer from the tables."
    display_and_speak(please_select_primer)
    while True:
        try:
            forward_index = int(input("Enter the index of the forward primer: ")) - 1
            reverse_index = int(input("Enter the index of the reverse primer: ")) - 1
        except ValueError:
            invalid_input = "Invalid input. Please enter a number."
            display_and_speak(invalid_input)
            continue
        else:
            try:
                forward_primer = forward_primers[forward_index]
                reverse_primer = reverse_primers[reverse_index]
                if reverse_complement(forward_primer[-3:]) == (reverse_primer[:3]):
                    unsuitable_basic = "The selected primers are not suitable as a primer pair."
                    display_and_speak(unsuitable_basic)
                elif abs(forward_temps[forward_index] - reverse_temps[reverse_index]) > 5:
                    unsuitable_temps = "The selected primers are not suitable as a primer pair, melting temperatures are too different. >5°C"
                    display_and_speak(unsuitable_temps)
                else:
                    suitable = "The selected primers are suitable as a primer pair."
                    display_and_speak(suitable)
                    break
            except IndexError:
                invalid_index = "Invalid index. Please enter a valid index."
                display_and_speak(invalid_index)
                continue

        primer_similarity = levenshtein_distance_algo(forward_primer, reverse_primer)
        print(f"The similarity of the complementary primer sequences is {primer_similarity}%")


def levenshtein_distance_algo(forward, reverse):
    lev_distance = distance(reverse_complement(forward), reverse)
    similarity = (1 - (lev_distance / max(len(forward), len(reverse)))) * 100
    return similarity


def save_all_files(gene_accession, snp_position, forward_tables, reverse_tables, sequence_data):
    try:
        with open(f"{gene_accession}_snp_{snp_position}_primers.txt", "w") as file:
            file.write(forward_tables + "\n" + reverse_tables)
        with open(f"{gene_accession}_snp_{snp_position}_sequence_data.txt", "w") as file:
            file.write(sequence_data)
    except IOError as e:
        saving_error = f"An error occurred while trying to save the files: {e}"
        check_permissions = "Please check if you have write permissions in the directory or if there is sufficient disk space."
        display_and_speak(saving_error, check_permissions)


if __name__ == "__main__":
    main()
