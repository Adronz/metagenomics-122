import os
import csv
import time



def fasta_to_string(file_path):
    with open(file_path, 'r') as file:
        # Skip header lines and concatenate sequence lines
        return ''.join(line.strip() for line in file if not line.startswith('>'))
    
    
def fasta_to_string_w_header(file_path):
    with open(file_path, 'r') as file:
        # Skip header lines and concatenate sequence lines
        return ''.join(line.strip() for line in file)

def txt_to_string(file_path):
    with open(file_path, 'r') as file:
        file_contents = file.read()
        return file_contents
        
def fasta_to_dict(file_path):
    fasta_dict = {}
    current_header = ""

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):  # Header line
                current_header = line[1:]  # Remove the '>' character
                fasta_dict[current_header] = ""
            else:  # Sequence line
                fasta_dict[current_header] += line

    return fasta_dict


#########################
#The Real Data For 1D
#########################

# Specify the path to your folder
folder_path = '/Users/AdrianHanson/Downloads/project1d (1)'

# List all files in the folder
file_names = os.listdir(folder_path)

# Filter out directories, keeping only files
file_paths = [os.path.join(folder_path, file_name) for file_name in file_names if os.path.isfile(os.path.join(folder_path, file_name))]

gametime_reads = None
for i in file_paths:
    if i == '/Users/AdrianHanson/Downloads/project1d (1)/project1d_reads.fasta':
        gametime_reads = i
        file_paths.remove(i)
        break


gametime_reads = fasta_to_dict(gametime_reads)
# print(gametime_reads)
#sort all genome names 

def extract_genome_number(file_path):
    # Split by '_', take the last part, and then split by '.' to get the number
    return int(file_path.split('_')[-1].split('.')[0])

# Sort the list by genome number
sorted_paths = sorted(file_paths, key=extract_genome_number)

#make gametime paths
gametime_genomes = []
for path in sorted_paths:
    new_genome = fasta_to_string(path)
    gametime_genomes.append(new_genome)
    
# for file in file_paths != '/Users/AdrianHanson/Downloads/project1c/project1c_reads.fasta':

#############################
#############################
start_time = time.time()

# Custom sorting function
def extract_number(s):
    # Extracting the number part and converting it to integer
    return int(s.split()[0][2:])


    #helper functions to make building this out easier
def read_minimizer(w:int, k:int, single_read:str,):
    read_hash = minimizer_sketch(w, k, single_read)
    return read_hash

def genome_minimizers(w:int,k:int,genome:str):
    genome_hash = minimizer_sketch(w, k, genome)
    return genome_hash

def hamming_distance(read:str, genome:str, aligned_loc:int):
    dist = 0
    for i in range(len(read)):
        if aligned_loc + i >= len(genome):
            break
        if read[i] != genome[aligned_loc + i]:
            dist +=1
    return dist

def minimizer_sketch(w:int,k:int,seq:str):
    N = len(seq)
    min_hash = {}
    last_min = 'G'*100
    last_loc = -1
    for i in range(N-k+1): #go through the genome or read
        temp_min = seq[i:i+k]
        temp_loc = i
        window = min(i + w - k, N)
        #window = i + w - k #location in seq + window length - kmer length so it stops k bp before the end 
        for j in range(i,window): #find the minimum kmer inside the window
            
            if last_loc < i: #if the old minimizer is outside of the window, automatically update to a minimizer 
                    last_min = seq[j:j+k]
                    last_loc = j
                    temp_min = seq[j:j+k]
                    temp_loc = j

            #everything else should be inside the window 
            if seq[j:j+k] < temp_min:
                temp_min = seq[j:j+k]
                temp_loc = j
        
        #now compare the minimum of the new window to the old minimizer 
        if temp_min < last_min: #if new minimizer is found
             last_min = temp_min 
             last_loc = temp_loc

             if last_min not in min_hash:
                min_hash[last_min] = [last_loc]
             elif last_min in min_hash:
                min_hash[last_min].append(last_loc)

    return min_hash




#finds the best match for a read in the reference, does not account for indels
def is_in_genome(ref_genome:str, full_read:str, genome_hash:{}, read_hash:{}, dist_threshold):
    min_dist = 10000
    for key in read_hash: #this iterates through the hash of each read 
        read_list = read_hash[key]
        read_length = len(key)
        #account for duplicate minimizers in read
        for read_min_loc in read_list:
            if key in genome_hash: #if the minimizer is also a genome minimizer 
                genome_loc_list = genome_hash[key] #get all possible locations
                for i in genome_loc_list: #get hamming distance for each possible match
                    aligned_loc = i - read_min_loc
                    temp_dist = hamming_distance(full_read,ref_genome, aligned_loc)
                    if temp_dist < min_dist:
                        min_dist = temp_dist

    if min_dist < dist_threshold:
        return True
    else:
        return False


def first_pass(w:int,k:int, genomes:[], reads:{}, tolerance:int):
#Idea: if all minimizers, or minimizers within a tolerance match the genome, add one to that genomes score.

    genome_scores = dict()

    read_hashes = [] #lets not remake anything every time thats stupid. 40,000,000 hopefully wont go to disk
    for read in reads: 
        read_hash = read_minimizer(w,k,reads[read])
        read_hashes.append(read_hash)

    for i in range(len(genomes)): #go through genomes one by one, store all read hashes
        num_minimizers = None
        genome_mins = minimizer_sketch(w,k,genomes[i]) 
        iter_str = str(i) #convert iterator to a string so it can be used as a key
        for read_hash in read_hashes:
            match_counter = 0
            num_minimizers = len(read_hash)
            for minimizer in read_hash:
                if minimizer in genome_mins:
                    match_counter += 1
        
        if match_counter >= num_minimizers - tolerance:
            if iter_str not in genome_scores:
                genome_scores[iter_str] = 1
            else:
                genome_scores[iter_str] += 1

    end_time = time.time()
    runtime = end_time - start_time
    print(runtime/60, "minutes")
    return genome_scores


def first_pass(w:int,k:int, genomes:[], reads:{}, tolerance:int):
#Idea: if all minimizers, or minimizers within a tolerance match the genome, add one to that genomes score.

    genome_scores = dict()

    read_hashes = [] #lets not remake anything every time thats stupid. 40,000,000 hopefully wont go to disk
    for read in reads: 
        read_hash = read_minimizer(w,k,reads[read])
        read_hashes.append(read_hash)

    for i in range(len(genomes)): #go through genomes one by one, store all read hashes
        num_minimizers = None
        genome_mins = minimizer_sketch(w,k,genomes[i]) 
        iter_str = str(i) #convert iterator to a string so it can be used as a key
        # print(genome_mins)

        for read_hash in read_hashes:
            # print(read_hash)
            match_counter = 0
            num_minimizers = len(read_hash)
            # print(num_minimizers, "minimizers")
            for minimizer in read_hash:
                if len(minimizer) < k:
                    continue
                # print(minimizer)
                if minimizer in genome_mins:
                    match_counter += 1
            # print("number of matches: ", match_counter)
          
            if match_counter >= num_minimizers - tolerance:
                if iter_str not in genome_scores:
                    genome_scores[iter_str] = 1
                else:
                    genome_scores[iter_str] += 1
        
    end_time = time.time()
    runtime = end_time - start_time
    print(runtime/60, "minutes")
    return genome_scores


w = 25
k = 10
threshold = 5

genome_tester = gametime_genomes[:100]




#reads_in_genomes_hash = all_together_now(w,k, genome_tester, gametime_reads, threshold)

genomes_present = first_pass(w, k, gametime_genomes, gametime_reads, 0)
print("there are ", len(genomes_present)," genomes present in the sample with tolerance 0")

file_name = '1d_genomes-0_dict.csv'

# Open the file for writing
with open(file_name, mode='w', newline='') as file:
    writer = csv.writer(file)
    
    # Write the header (optional)
    writer.writerow(['Key', 'Value'])
    
    # Write the dictionary items (keys and values)
    for key, value in genomes_present.items():
        writer.writerow([key, value])


genomes_present2 = first_pass(w, k, gametime_genomes, gametime_reads, 1)
print("there are ", len(genomes_present)," genomes present in the sample with tolerance 1")

file_name2 = '1d_genomes-1_dict.csv'
# Open the file for writing
with open(file_name2, mode='w', newline='') as file:
    writer = csv.writer(file)
    
    # Write the header (optional)
    writer.writerow(['Key', 'Value'])
    
    # Write the dictionary items (keys and values)
    for key, value in genomes_present2.items():
        writer.writerow([key, value])




# prediction_list = []
# for hash in reads_in_genomes_hash:
#     this_hash = reads_in_genomes_hash[hash]
#     for i in range(len(this_hash)):
#         prediction_list.append(this_hash[i])


# Writing to the CSV file
# with open('predictions.csv', 'w', newline='') as csv_file:
#     # Step 4: Using csv.writer to write the list to the CSV file
#     writer = csv.writer(csv_file)
#     writer.writerow(prediction_list) # Use writerow for single list


# with open("predictions.csv", "w", newline='') as csv_file:
#     csv_writer = csv.writer(csv_file)
#     csv_writer.writerow(["predictions"])

#     for snp in prediction_list:
#      csv_writer.writerow([snp])




# for key in reads_in_genomes_hash:
#     this_shit = reads_in_genomes_hash[key]
#     for read in this_shit:
#         print(read)

# matches = 0
# denominator = len(genome_all_soln)
# for hash in reads_in_genomes_hash:
#     this_hash = reads_in_genomes_hash[hash]
#     for read in this_hash:
#         # print(read)
#         if read in genome_all_soln:
#              matches += 1

# print("there are ", matches,"matches" )
# print((matches/denominator)*100, "% correct")



# genome_0 = fasta_to_string(genome_0)
# genome_1 = fasta_to_string(genome_1)
# genome_2 = fasta_to_string(genome_2)
# genome_3 = fasta_to_string(genome_3)
# genome_4 = fasta_to_string(genome_4)
# genome_5 = fasta_to_string(genome_5)
# genome_6 = fasta_to_string(genome_6)
# genome_7 = fasta_to_string(genome_7)
# genome_8 = fasta_to_string(genome_8)
# genome_9 = fasta_to_string(genome_9)

# ordered_genomes12 = [genome_0 ,genome_1, genome_2, genome_3, genome_4, genome_5, genome_6, genome_7, genome_8, genome_9]
# # remainder = [ genome_1, genome_2, genome_3, genome_4, genome_5, genome_6, genome_7, genome_8, genome_9]

# sample_reads = fasta_to_dict(sample_reads)

# genome_with_source = fasta_to_string_w_header(genome_with_source)
# print(genome_with_source[0:1000])
#all of these should return true
# genome_1a = fasta_to_string(project_1a_data)
# reads_1a = fasta_to_dict(project_1a_reads)
# solutions = txt_to_string(solutions)
# solutions = solutions.split('\n')

# for i in range(0, 100):
#      print(solutions[i])


# Take only the first part (before the tab)
    
################################
#TESTING PURPOSES
################################

#This checks if its in the whole genome which is like, duh. It needs to check if the reads are in genome 2 vs genome 1. 
#maybe run each separately /

#see how well my program did with genome 0 alone
# genome_all_soln = []
# for solution in solutions:
#         if "Genome_Number_0" in solution:
#             genome_all_soln.append(str(solution))

# for i in range(len(genome_all_soln)):
#     # Split each element by the tab character ('\t')
#     parts = genome_all_soln[i].split('\t')

#     # Take only the first part (before the tab) and update the element in the list
#     genome_all_soln[i] = parts[0]

###############################################################
###############################################################
