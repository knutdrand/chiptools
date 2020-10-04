import click
import matplotlib.pyplot as plt
import seaborn as sns
from .pandas_parsers import read_bedgraph, read_bedfile
from .aggregateplot import * #TSSPlot, VPlot, HeatPlot, AveragePlot
plot_types = {"v": VPlot, "average": AveragePlot, "heat": HeatPlot, "tss": TSSPlot}
sns.set_theme()

@click.command()
@click.argument("plot_type", type=str)
@click.argument("bedgraph", type=click.File("r"))
@click.argument("bedfile", type=click.File("r"))
@click.option("-o", "--out_im", "out_im", type=click.File("wb"))
@click.option("-od", "--out_data", "out_data", type=click.File("wb"))
@click.option("-w", "--width", "figure_width", default=2000)
def main(plot_type, bedgraph, bedfile, out_im, out_data, figure_width): #, figure_width, region_size, aspect_ratio):
    assert plot_type in plot_types
    bedgraphs = read_bedgraph(bedgraph)
    regions = read_bedfile(bedfile)
    f = plot_types[plot_type](figure_width=figure_width)
    fig = f(bedgraphs, regions)
    if "x" in fig:
        p = sns.lineplot(data=fig, x="x", y="y")
    else:
        p = sns.heatmap(fig)
    p.set_xlabel(f.xlabel)
    p.set_ylabel(f.ylabel)
    p.set_title(f"{plot_type}-plot")
    if out_im is not None:
        plt.savefig(out_im)
    if out_data is not None:
        fig.to_pickle(out_data)
    if (out_data is None) and (out_im is None):
        plt.show()
