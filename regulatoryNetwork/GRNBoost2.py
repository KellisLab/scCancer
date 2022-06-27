# import python modules
import sys
import pandas as pd
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2, genie3
from distributed import LocalCluster, Client

if __name__ == '__main__':
    wd = sys.argv[1]
    # load the data
    net1_tf_path = wd + 'int/1.1_inputTFs.txt'
    net1_ex_path = wd + 'int/1.1_exprMatrix_filtered_t.txt'

    ex_matrix = pd.read_csv(net1_ex_path, sep='\t')
    tf_names = load_tf_names(net1_tf_path)

    local_cluster = LocalCluster(n_workers=16, threads_per_worker=1)

    custom_client = Client(local_cluster)

    # infer the gene regulatory network
    network = grnboost2(expression_data=ex_matrix, tf_names=tf_names, client_or_address=custom_client)

    net1_out_path = wd + 'int/1.1_network.tsv'

    network.to_csv(net1_out_path, sep='\t', header=False, index=False)

    custom_client.close()
    local_cluster.close()
