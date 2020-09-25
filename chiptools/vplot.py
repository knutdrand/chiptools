import numpy as np

def vplot(bedgraph, regions, size_x, size_y, max_size):
    sizes = regions.ends-regions.starts
    diffs = np.zeros((size_y, size_x))
    mids = (regions.ends+regions.starts)//2
    signals = bedgraph.get_slices(mids-max_size//2, mids+max_size//2, regions.directions)
    Ns = np.zeros(size_y)
    for size, signal in zip(sizes, signals):
        if size>max_size:
            continue
        assert np.all(signal._values >= 0)
        graph_diffs = signal.scale_x(size_x).to_graph_diffs()
        graph_diffs.assert_positive()
        row = min(int(size/max_size*size_y), diffs.shape[0]-1)
        assert np.all(np.cumsum(diffs[row])>=0), np.where(np.cumsum(diffs[row])<0)
        graph_diffs.update_dense_array(diffs[row])
        assert np.all(np.cumsum(diffs[row])>=0), (graph_diffs, np.cumsum(diffs[row]))
        Ns[row]+=1
    tot = np.cumsum(diffs, axis=1)
    assert np.all(tot>=0)
    return tot, Ns

def get_heatplot(regions, bedgraphs):
    sizes = np.zeros(sum(r.starts.size for r in regions.values()))
    indices = {}
    cur_index = 0
    for chrom, r in regions.items():
        n = r.starts.size
        sizes[cur_index:cur_index+n] = r.ends-r.starts
        indices[chrom]=(cur_index, cur_index+n)
        cur_index+=n

    # sizes = np.concatenate([regions[chrom].ends-r[chrom].starts for chrom in bedgraphs
    #                         if chrom in regions and chrom!="chrM"])
    args_tmp = np.argsort(sizes)
    args = np.empty_like(args_tmp)
    args[args_tmp] = np.arange(len(args))
    N = args.size
    size_x = 10000
    size_y = 20000
    max_size=50000
    signal = np.zeros((size_y, size_x))
    Ns = np.zeros(size_y)
    for chrom, bedgraph in bedgraphs:
        if chrom not in regions or chrom=="chrM":
            continue
        assert "alt" not in chrom, chrom
        print("Reading", chrom)
        chrom_regions = regions[chrom]

        cur_args=args[indices[chrom][0]:indices[chrom][1]]
        # args_index+=chrom_regions.starts.size
        s, Ns_tmp = heatplot_per_chrom(bedgraph, chrom_regions, size_x, size_y, max_size, cur_args, N)
        signal += s
        Ns += Ns_tmp
    singal = signal/(np.maximum(Ns, 1)[:, None])
    return signal


def heatplot_per_chrom(bedgraph, regions, size_x, size_y, max_size, args, N):
    sizes = regions.ends-regions.starts
    indices = np.argsort(sizes)
    diffs = np.zeros((size_y, size_x))
    mids = (regions.ends+regions.starts)//2
    signals = bedgraph.get_slices(mids-max_size//2, mids+max_size//2, regions.directions)
    Ns = np.zeros(size_y)
    for arg, signal in zip(args, signals):
        graph_diffs = signal.to_graph_diffs().scale_x(size_x)
        graph_diffs.update_dense_array(diffs[int(arg/N*size_y)])
        Ns[int(arg/N*size_y)]+=1
    return np.cumsum(diffs, axis=1), Ns
