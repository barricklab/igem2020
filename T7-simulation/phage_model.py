from Bio import SeqIO
import pinetree as pt
import sys
import os
import datetime
import argparse
import multiprocessing
import random

CELL_VOLUME = 1.1e-15
PHI10_BIND = 1.82e7  # Binding constant for phi10


IGNORE_REGULATORY = ["E. coli promoter E[6]",
                     "T7 promoter phiOR",
                     "T7 promoter phiOL",
                     "E. coli promoter A0 (leftward)"]

### Some important proteins
#
# gp1 is T7 RNA Polymerase
# gp3.5 is T7 lysozyme
# gp0.7 is a protein kinase that phosphorylates E. coli RNA polymerase

#New 'CDS' feature version
IGNORE_CDS = ["gp10B",
              "gp5.5-5.7",
              "gp4.1",
              "gp4B",
              "gp0.6A",
              "gp0.6B",
              "gp0.5",
              "gp0.4"]


class Logger:
    '''Sends pretty colors to the console and also logs console to file'''
    def __init__(self, log_output = "", verbose = False):
        self.colors = {'normal': "\u001b[0m",
                  'warn': '\u001b[31m'}
        self.verbose = verbose
        self.log_output = log_output

        if self.log_output:  # Gotta make sure this exists
            self.log_output = self.log_output.replace("\\", "/")
            out_dir = "/".join(self.log_output.split("/")[:-1])
            if (out_dir != '') and (not os.path.exists(out_dir)):
                os.makedirs(out_dir)
            if self.log_output[-1] == "/" or self.log_output[-1] == ".":
                self.log_output = f"{self.log_output}pinetree.log"
            with open(self.log_output, "w") as _:  # clears the file
                pass

    def _send_to_log(self, text):
        if not self.log_output:
            return
        with open(self.log_output, 'a') as file:
            file.write(text + '\n')


    def normal(self,text):
        print(f"{self.colors['normal']}{text}{self.colors['normal']}")
        self._send_to_log(f"[NORMAL] {text}")
    def warn(self, text):
        print(f"{self.colors['warn']}Warning: {text}{self.colors['normal']}")
        self._send_to_log(f"[WARNING] {text}")
    def log(self, text):
        if self.verbose:
            print(f"{self.colors['normal']}{text}{self.colors['normal']}")
        self._send_to_log(f"[LOG] {text}")


# Optimal E. coli codons
OPT_CODONS_E_COLI = {'A': ['GCT'],
                     'R': ['CGT', 'CGC'],
                     'N': ['AAC'],
                     'D': ['GAC'],
                     'C': ['TGC'],
                     'Q': ['CAG'],
                     'E': ['GAA'],
                     'G': ['GGT', 'GGC'],
                     'H': ['CAC'],
                     'I': ['ATC'],
                     'L': ['CTG'],
                     'F': ['TTC'],
                     'P': ['CCG'],
                     'S': ['TCT', 'TCC'],
                     'T': ['ACT', 'ACC'],
                     'Y': ['TAC'],
                     'V': ['GTT', 'GTA']}


def get_promoter_interactions(name):
    '''
    Calculate promoter binding strengths. The relative strengths defined here
    come from 2012 Covert, et al paper.
    '''
    if name in IGNORE_REGULATORY:
        return
    ecoli_strong = ["E. coli promoter A1",
                    "E. coli promoter A2",
                    "E. coli promoter A3"]
    ecoli_weak = ["E. coli B promoter",
                  "E. coli C promoter"]
    phi1_3 = ["T7 promoter phi1.1A",
              "T7 promoter phi1.1B",
              "T7 promoter phi1.3",
              "T7 promoter phi1.5",
              "T7 promoter phi1.6"]
    phi3_8 = ["T7 promoter phi2.5",
              "T7 promoter phi3.8",
              "T7 promoter phi4c",
              "T7 promoter phi4.3",
              "T7 promoter phi4.7"]
    phi6_5 = ["T7 promoter phi6.5"]
    phi9 = ["T7 promoter phi9"]
    phi10 = ["T7 promoter phi10"]
    phi13 = ["T7 promoter phi13",
             "T7 promoter phi17"]

    if name in ecoli_strong:
        return {'ecolipol': 10e4,
                'ecolipol-p': 3e4}
    elif name in ecoli_weak:
        return {'ecolipol': 1e4,
                'ecolipol-p': 0.3e4}
    elif name in phi1_3:
        return {'gp1': PHI10_BIND * 0.01,
                'gp3.5': PHI10_BIND * 0.01 * 0.5}
    elif name in phi3_8:
        return {'gp1': PHI10_BIND * 0.01,
                'gp3.5': PHI10_BIND * 0.01 * 0.5}
    elif name in phi6_5:
        return {'gp1': PHI10_BIND * 0.05,
                'gp3.5': PHI10_BIND * 0.05}
    elif name in phi9:
        return {'gp1': PHI10_BIND * 0.2,
                'gp3.5': PHI10_BIND * 0.2}
    elif name in phi10:
        return {'gp1': PHI10_BIND,
                'gp3.5': PHI10_BIND}
    elif name in phi13:
        return {'gp1': PHI10_BIND * 0.1,
                'gp3.5': PHI10_BIND * 0.1}
    else:
        raise ValueError(
            "Promoter strength for {0} not assigned.".format(name))


def get_terminator_interactions(name):
    '''
    Get terminator efficiencies.
    '''
    if name == "E. coli transcription terminator TE":
        return {'ecolipol': 1.0,
                'ecolipol-p': 1.0,
                'gp1': 0.0,
                'gp1+gp3.5': 0.0}
    elif name == "T7 transcription terminator Tphi":
        return {'gp1': 0.85,
                'gp1+gp3.5': 0.85}
    else:
        return {'name': 0.0}


def compute_cds_weights(record, feature, factor, weights):
    # Grab the gene name
    nuc_seq = feature.location.extract(record).seq
    aa_seq = feature.qualifiers["translation"][0]
    weight_sum = 0
    for index, nuc in enumerate(nuc_seq):
        aa_index = int(index / 3)
        codon_start = aa_index * 3
        codon = nuc_seq[codon_start:codon_start + 3]
        genome_index = feature.location.start + index
        if aa_index < len(aa_seq):
            if aa_seq[aa_index] in OPT_CODONS_E_COLI:
                if codon in OPT_CODONS_E_COLI[aa_seq[aa_index]]:
                    weights[genome_index] = factor
                    weight_sum += factor
                else:
                    weights[genome_index] = 1
                    weight_sum += 1
    return weights


def normalize_weights(weights):
    # Average over all CDSs, which will have non-zero weights
    non_zero = sum(1 if i != 0 else 0.0 for i in weights)
    mean_weight = sum(weights) / non_zero
    norm_weights = [i / mean_weight for i in weights]
    # Replace non-CDS weights with 1
    norm_weights = [1 if i == 0 else i for i in norm_weights]
    return norm_weights


def phage_model(input, output=None, time=1500, verbose=True, seed=None):

    sim = pt.Model(cell_volume=CELL_VOLUME)

    record = SeqIO.read(input, "genbank")
    genome_length = len(record.seq)
    phage = pt.Genome(name="phage", length=genome_length)

    if not output:
        output = ".".join(input.split(".")[:-1])
    # Make the directory for output if it doesnt exist
    output_dir = output.replace("\\", "/")
    output_dir = "/".join(output.split("/")[:-1])
    if output_dir != '' and not os.path.exists(output_dir):
      os.makedirs(output_dir)


    # Log relevant information
    if output[-1] == "/" or output[-1] == ".":
        log_output = f"{output}pinetree.log"
    else:
        log_output = f"{output}.log"
    logger = Logger(log_output=f"{log_output}", verbose=verbose)
    start_time = datetime.datetime.utcnow()
    logger.normal("[Pinetree] Pinetree T7 Genome Simulation")
    logger.normal("barricklab/igem2020 Fork")
    # Try and find a git repo and log its last commit
    if os.path.exists(".git/refs/heads/master"):
        git_master_path =  ".git/refs/heads/master"
    elif os.path.exists("../.git/refs/heads/master"):
        git_master_path = "../.git/refs/heads/master"
    else:
        git_master_path = ""
    if git_master_path:
        with open(git_master_path, 'r') as file:
            commit_hash = file.readline().strip()
            logger.normal(f"Last commit: {commit_hash}")
    logger.normal(f"Script and simulation started at {start_time} UTC")

    # --- Feature Acquisition and validation  VVV

    feature_dict = {}
    for i, feature in enumerate(record.features):  # Accuasition
        start = feature.location.start.position + 1
        stop = feature.location.end.position
        name = ''
        feature_type = ''
        interactions = None
        source_feature = feature
        if 'name' in feature.qualifiers:
            name = feature.qualifiers["name"][0]
        elif "note" in feature.qualifiers:
            name = feature.qualifiers["note"][0]
        if feature.type == "regulatory":
            if "promoter" in feature.qualifiers["regulatory_class"]:
                length = stop - start
                if length < 35:
                    start = start - 35
                interactions = get_promoter_interactions(name)
                feature_type = "promoter"
            if "terminator" in feature.qualifiers["regulatory_class"]:
                interactions = get_terminator_interactions(name)
                feature_type = "terminator"
        elif feature.type == "CDS":
            feature_type = "cds"
        else:
            feature_type = None

        if feature_type:
            feature_dict[i] = {"start": start,
                               "stop": stop,
                               "name": name,
                               "type": feature_type,
                               "interactions": interactions,
                               "skip": False,
                               "source_feature": source_feature
                               }

    # TODO: Add more feature validation
    for feature in feature_dict.items():  # Validation
        feature = feature[1]
        if feature['name'] in IGNORE_REGULATORY or feature['name'] in IGNORE_CDS:
            feature['skip'] = True
            logger.log(f"Ignored feature {feature['name']} ({feature['start']} - {feature['stop']})")
            continue
        if 'pseudo' in feature['source_feature'].qualifiers.keys():
            logger.warn(f"Found {feature['name']} with flag 'pseudo'. Ignoring.")
            feature['skip'] = True
            continue
        if feature['stop'] - feature['start'] < 50 and feature['type'] in ['gene', 'cds']:
            logger.warn(f"Found {feature['type'], feature['name']} that is tiny! ({feature['start']} - {feature['stop']})")

    # -- Feature Acquisition Validation  ^^^
    # -- Add Features to Sim  VVV
    for feature in feature_dict.items():
        feature = feature[1]
        weights = [0.0] * len(record.seq)  # TODO: Maybe this is a bug and should be outside loop?
        if feature['skip']:
            continue
        elif feature['type'] == "promoter":
            phage.add_promoter(feature['name'], feature['start'], feature['stop'], feature['interactions'])
            logger.log(f"Added promoter feature: {feature['name']}, Start: {feature['start']}, Stop: {feature['stop']}")
        elif feature['type'] == "terminator":
            phage.add_terminator(feature['name'], feature['start'], feature['stop'], feature['interactions'])
            logger.log(f"Added terminator feature: {feature['name']}, Start: {feature['start']}, Stop: {feature['stop']}")
        elif feature['type'] == "cds":
            phage.add_gene(name=feature['name'], start=feature['start'], stop=feature['stop'],
                           rbs_start=feature['start'] - 30, rbs_stop=feature['start'], rbs_strength=1e7)
            weights = compute_cds_weights(record, feature['source_feature'], 1.0, weights)
            logger.log(f"Added CDS feature: {feature['name']}, Start: {feature['start']}, Stop: {feature['stop']}")
        else:
            continue
    # -- Add Featues to Sim  ^^^



    logger.normal("Registered genome features")

    mask_interactions = ["gp1", "gp1+gp3.5",
                         "ecolipol", "ecolipol-p", "ecolipol-2", "ecolipol-2-p"]
    phage.add_mask(500, mask_interactions)

    norm_weights = normalize_weights(weights)
    phage.add_weights(norm_weights)

    logger.normal("Implemented masks and weighting")

    sim.register_genome(phage)

    sim.add_polymerase("gp1", 35, 230, 0)
    sim.add_polymerase("gp1+gp3.5", 35, 230, 0)
    sim.add_polymerase("ecolipol", 35, 45, 0)
    sim.add_polymerase("ecolipol-p", 35, 45, 0)
    sim.add_polymerase("ecolipol-2", 35, 45, 0)
    sim.add_polymerase("ecolipol-2-p", 35, 45, 0)

    sim.add_ribosome(30, 30, 0)

    sim.add_species("bound_ribosome", 10000)

    sim.add_species("bound_ecolipol", 1800)
    sim.add_species("bound_ecolipol_p", 0)
    sim.add_species("ecoli_genome", 0)
    sim.add_species("ecoli_transcript", 0)

    sim.add_reaction(1e6, ["ecoli_transcript", "__ribosome"], [
                     "bound_ribosome"])

    sim.add_reaction(0.04, ["bound_ribosome"], [
                     "__ribosome", "ecoli_transcript"])

    sim.add_reaction(0.001925, ["ecoli_transcript"], ["degraded_transcript"])

    sim.add_reaction(1e7, ["ecolipol", "ecoli_genome"], ["bound_ecolipol"])

    sim.add_reaction(
        0.3e7, ["ecolipol-p", "ecoli_genome"], ["bound_ecolipol_p"])

    sim.add_reaction(0.04, ["bound_ecolipol"], [
                     "ecolipol", "ecoli_genome", "ecoli_transcript"])

    sim.add_reaction(0.04, ["bound_ecolipol_p"], [
                     "ecolipol-p", "ecoli_genome", "ecoli_transcript"])

    sim.add_reaction(3.8e7, ["gp0.7", "ecolipol"],
                     ["ecolipol-p", "gp0.7"])

    sim.add_reaction(3.8e7, ["gp0.7", "ecolipol+gp2"],
                     ["ecolipol+gp2-p", "gp0.7"])

    sim.add_reaction(3.8e7, ["gp2", "ecolipol"], ["ecolipol+gp2"])

    sim.add_reaction(3.8e7, ["gp2", "ecolipol-p"], ["ecolipol+gp2-p"])

    sim.add_reaction(1.1, ["ecolipol+gp2-p"], ["gp2", "ecolipol-p"])

    sim.add_reaction(1.1, ["ecolipol+gp2"], ["gp2", "ecolipol"])

    sim.add_reaction(3.8e9, ["gp3.5", "gp1"], ["gp1+gp3.5"])

    sim.add_reaction(3.5, ["gp1+gp3.5"], ["gp3.5", "gp1"])

    logger.normal("Registered reactions")

    logger.normal("Running simulation")

    if not seed:
        seed = random.randint(0, 10000)
    sim.seed(seed)
    logger.normal(f"Random seed was set to {seed}")

    if output[-1] == "/" or output[-1] == ".":
        sim_output = f"{output}phage.counts.tsv"
    else:
        sim_output = f"{output}.counts.tsv"

    # -- Running the actual sim. VVVV
    # Note: Multiprocessing necessary for working keyboard interrupts.
    try:
        sim_process = multiprocessing.Process(target=sim.simulate, kwargs={'time_limit': time,
                                                                           'time_step': 5,
                                                                           'output': sim_output
                                                                           })
        sim_process.start()
        sim_process.join()
    except KeyboardInterrupt:
        sim_process.terminate()
        with open(sim_output, 'r') as outfile:
            for line in outfile:
                pass
            last_file_line = line
        interrupt_time = str(last_file_line).split("\t")[0]
        logger.warn(f'Received keyboard interruption. Simulation reached time {interrupt_time}')
        finish_time = datetime.datetime.utcnow()
        run_time = (finish_time - start_time).total_seconds()
        logger.normal(f"Simulation interrupted after {run_time / 60} minutes.")
        exit(0)
    except TypeError:
        logger.warn("There was a problem with multiprocessing. Keyboard Interrupts won't work, but the simulation should still run.")
        sim.simulate(time_limit= time, time_step=5, output=sim_output)

    finish_time = datetime.datetime.utcnow()
    run_time = (finish_time-start_time).total_seconds()

    logger.normal(f"Simulation completed in {run_time/60} minutes.")
    # -- Running the actual sim. ^^^^

if __name__ == "__main__":
    arguments = sys.argv
    # For hard coding variables if you feel like it:
    input_genome = None  # ex. (resources/T7_genome.gb)
    output_path = None   # ex. [output | output/]

    # Otherwise it will take from command line

    parser = argparse.ArgumentParser(description='Perform simulation of T7 phage gene expression')
    parser.add_argument('-i',
                       action='store',
                       metavar="input.gb",
                       dest='i',
                       required=True,
                       type=str,
                       help="input file in fasta format (REQUIRED)")

    parser.add_argument('-o',
                        action='store',
                        metavar="output-prefix",
                        dest='o',
                        required=False,
                        type=str,
                        help="prefix of *.counts.tsv and *.log output files (REQUIRED)")
    parser.add_argument('-t',
                        action='store',
                        metavar="seconds",
                        dest='t',
                        required=False,
                        type=int,
                        help="Duration of simulation in seconds")
    parser.add_argument('-s',
                        action='store',
                        metavar="seed",
                        dest='s',
                        required=False,
                        type=int,
                        help="Randomness seed")
    options = parser.parse_args()
    try:
        output_path = options.o
    except AttributeError:
        pass
    if options.i:
        input_genome = options.i
    try:
        time = options.t
    except AttributeError:
        time = 1500
    if not time:
        time = 1500

    try:
        seed = options.s
    except AttributeError:
        seed = None
    if not time:
        seed = None

    if not output_path:
        output_path = ".".join(input_genome.split(".")[:-1])

    if not os.path.exists(input_genome):
        print(f"Could not find file {input_genome}")
        exit(1)

    phage_model(input_genome, output_path, time, verbose=False, seed=seed)
