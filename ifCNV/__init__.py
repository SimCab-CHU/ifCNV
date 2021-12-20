#
# Copyright (C) 2021
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

__author__ = 'Simon Cabello'
__copyright__ = 'Copyright (C) 2021'
__license__ = 'GNU General Public License'
__version__ = '0.2.1'
__email__ = 's-cabelloaguilar@chu-montpellier.fr'
__status__ = 'prod'


import os
import re
import numpy as np
import pandas as pd
import subprocess
from sklearn.ensemble import IsolationForest
import plotly
import plotly.graph_objs as go
from pybedtools import BedTool


###########################################
#               FUNCTIONS                 #
###########################################


def generateReport(final, output_dir, reads, ressources, CNVneg):
    if not os.path.isdir(output_dir):
        cmd = ["mkdir", output_dir]
        subprocess.check_output(cmd)

    cmd = ["cp", "-r", ressources, output_dir]
    subprocess.check_output(cmd)

    output_report = f"{output_dir}/run.html"
    html = getTemplate()

    for i in final.index:
        fff = final.loc[i]

        run = fff['Run']
        sample = fff['Sample name']
        region = fff['Region']

        generateGraph(sample, output_dir, reads, region, CNVneg)

        ratio = str(round(fff['Reads ratio'], 2))
        score = str(round(fff['Score'], 2))

        html = html.replace('<!--ADD RUN-->', '<!--ADD RUN-->\n\t\t\t<!--ADD RUN-->')

        torep = (
            f'<tr url="{sample}_{region}.html" style="cursor: pointer;">'
            f'<td>{run}</td>'
            f'<td>{sample}</td>'
            f'<td>{region}</td>'
            f'<td>{ratio}</td>'
            f'<td>{score}</td></tr>'
        )
        html = html.replace('<!--ADD RUN-->', torep, 1)
        with open(output_report, "w") as FH_out:
            FH_out.write(html)


def generateGraph(sample, output_dir, reads, region, CNVneg):
    dfneg = np.mean(reads[CNVneg], axis=1)
    dfpos = reads[sample]

    df = np.log2(np.array(dfpos) / np.array(dfneg))

    genes = [i.split("_")[0] for i in reads.index]

    col = []
    for i in reads.index:
        if i.startswith(f"{region}_"):
            col.append("red")
        else:
            col.append("black")

    data = [go.Scatter(
        x=genes,
        y=df,
        mode='markers',
        marker=dict(color=col, size=10)
    )]
    layout = go.Layout(
        yaxis=dict(title='Nomalized reads ratio'),
        plot_bgcolor='rgb(240,240,240)'
    )
    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(
        fig,
        filename=f"{output_dir}/{sample}_{region}.html",
        auto_open=False
    )


def getTemplate():
    html = """
<!DOCTYPE html>
<html>
  <head>
    <title>ifCNV</title>
    <meta charset="UTF-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <link rel="stylesheet" href="ressources/lib///w3.css" />
    <link rel="stylesheet" href="ressources/fontawesome-free-5.3.1-web/css/all.css" />
    <link rel="stylesheet" href="ressources/lib///datatables/datatables.min.css" />
    <link rel="stylesheet" href="ressources/lib///interface.css" />
  </head>
  <body>
    <!-- Top container -->
    <div class="w3-bar w3-top w3-blue-grey w3-large" style="margin-left:10px;">
      <img src="ressources/ifCNV.png" width="170" height="80">
      <span style="color:white" class="w3-display-middle"><b><i><font size="+20" face="calibri">ifCNV Report</font></i></b></span>
    </div>
  <body>
    <div id="Summary" class="w3-container" style="margin-left:50px;margin-top:125px;">
      <h2>
        <b><center>Results</center></b>
      </h2>

    <div class="w3-main" style="margin-left:100px;margin-right:100px;">
        <div id="runs" class="w3-container" style="padding-top:10px;">
          <table id="T_runs" class="w3-table-all w3-hoverable w3-centered">
            <thead>
              <tr>
                <th>Run</th>
                <th>Sample Name</th>
                <th>Region</th>
                <th>Reads ratio</th>
                <th>Score</th>
              </tr>
            </thead>
            <tbody>
				<!--ADD RUN-->
            </tbody>
            <tfoot class="" style="">
              <tr class="" style="">
                <th>Run</th>
                <th>Sample Name</th>
                <th>Region</th>
                <th>Reads ratio</th>
                <th>Score</th>
              </tr>
            </tfoot>
          </table>
        </div>
    </div>
</body>

<script type="text/javascript" src="ressources/lib/datatables/datatables.min.js"></script>
  <script>
  $('body').on('mousedown', 'tr[url]', function(e){
      var click = e.which;
      var url = $(this).attr('url');
      if(url){
          if(click == 2 || (click == 1 && (e.ctrlKey || e.metaKey))){
              window.open(url, '_blank');
              window.focus();
          }
          else if(click == 1){
              window.location.href = url;
          }
          return true;
      }
  });

            $(document).ready(function() {
                $('#T_runs tfoot th').each( function () {
                    var title = $(this).text();
                    $(this).html( '<input type="text" placeholder="Search '+title+'" style="width:100%" />' );
                } );

                var table = $('#T_runs').DataTable({
                    paging:true,
                    order: [[1, 'desc']],
                });

                table.columns().every( function () {
                    var that = this;

                    $( 'input', this.footer() ).on( 'keyup change', function () {
                        if ( that.search() !== this.value ) {
                            that
                                .search( this.value )
                                .draw();
                        }
                    });
                });
            });
        </script>

</html>
    """
    return html


def replaceZeroes(data):
    min_nonzero = np.min(data[np.nonzero(data)])
    data[data == 0] = min_nonzero
    return data


def createReadsMatrix(pathToBam, bedFile, output=None, verbose=False):
    cmd = ["ls", pathToBam]
    res = subprocess.check_output(cmd)
    final = pd.DataFrame()

    bed = BedTool(bedFile)

    for i in res.decode('utf-8').split("\n"):
        if i.endswith(".bam"):
            sampleName = i[:-4]
            if verbose:
                print(f"Processing sample {sampleName}...")
            try:
                pathBam = pathToBam + "/" + i
                data = bed.multi_bam_coverage(bams=pathBam)
                df = data.to_dataframe()
                final[sampleName] = df[df.keys()[len(df.columns) - 1]]
                if verbose:
                    print(f"{sampleName} Done \n")
            except subprocess.CalledProcessError:
                print(f"{sampleName} : skipped")
                print("Hint: Index file (.bai) must be present in the folder \n")

    final.index = list(df["name"])

    if output is not None:
        if verbose:
            print("Reads matrix created !")
        final.to_csv(output, sep="\t")

    return final


def filterReads(reads, N, regtar=None, regsamp=None):
    col = reads.columns
    rows = reads.index
    if regtar is not None:
        reads = reads.filter(regex=regtar, axis=0)
    if regsamp is not None:
        reads = reads.filter(regex=f"^(?!{regsamp})")
    reads = reads.filter(regex="^(?!MSI)", axis=0)
    reads = reads.filter(regex="^(?!TN)")
    reads = reads.filter(regex="^(?!TP)")
    reads = reads.filter(regex="^(?!HD)")
    reads = reads.filter(regex="^(?!H2)")
    reads = reads.loc[reads.sum(axis=1) / len(reads.columns) > N, :]
    filtered_samples = col[~np.in1d(col, reads.columns)]
    filtered_targets = rows[~np.in1d(rows, reads.index)]
    return (reads, filtered_samples, filtered_targets)


def normalizeReads(reads):
    reads_norm = reads / reads.sum(axis=0)
    return reads_norm


def aberrantSamples(reads, conta='auto', verbose=True):
    if conta != 'auto' and conta != 'none':
        conta = float(conta)
    if conta == 'none':
        conta = 1 / reads.shape[1]

    tmp = np.percentile(reads, 99, axis=0) / np.mean(reads, axis=0)
    random_data = np.array(tmp).reshape(-1, 1)
    clf = IsolationForest(contamination=conta).fit(random_data)
    preds = clf.predict(random_data)
    res_amp = np.array(reads.columns)[preds == -1]

    tmp = np.percentile(reads, 1, axis=0) / np.mean(reads, axis=0)
    random_data = np.array(tmp).reshape(-1, 1)
    clf = IsolationForest(contamination=conta).fit(random_data)
    preds = clf.predict(random_data)
    res_del = np.array(reads.columns)[preds == -1]

    res = np.unique(np.concatenate((res_amp, res_del)))
    norm = np.array(reads.columns[~np.in1d(reads.columns, res)])

    if verbose:
        print("Aberrant sample(s): \n")
        for i in res:
            print(f"\t{i}\n")

        print("Normal sample(s): \n")
        for i in norm:
            print(f"\t{i}\n")

    return (res, norm)


def aberrantAmpliconsPerSample(name, reads_norm, CNVneg, conta=0.01):
    if conta != 'auto':
        if conta != 'none':
            conta = float(conta)
    if conta == 'none':
        conta = 1 / reads.shape[1]
    random_data = np.array(reads_norm[name])
    random_data = replaceZeroes(random_data)
    norm = np.array(np.mean(reads_norm[CNVneg], axis=1))
    norm = replaceZeroes(norm)
    df = np.log2(random_data / norm)
    clf = IsolationForest(contamination=conta).fit(df.reshape(-1, 1))
    preds = clf.predict(df.reshape(-1, 1))

    return np.array(reads_norm.index)[preds == -1]


def scoreAmplif(k, n, N):
    p = n / N
    x = np.log(1 / ((p**k) * (1 - p)**(n - k))) * (k / n)
    return x


def amplifEvalGene(reads_norm, region, CNVneg, sample):
    reads_m = reads_norm.filter(regex=f"^{region}", axis=0)
    pos = np.array(reads_m[sample])
    norm = np.array(np.mean(reads_m[CNVneg], axis=1))
    df = pos / norm
    val = np.mean(df)
    if val == np.inf:
        val = 100

    return val


def aberrantAmpliconsFinal(reads, reads_norm, CNVpos, CNVneg, scoreThreshold=10, conta=0.01, mode="fast", run="ifCNV", verbose=True):
    f = pd.DataFrame(columns=["Run", "Sample name", "Region", "Reads ratio", "Score"])

    if mode == "extensive":
        samples = [*CNVpos, *CNVneg]
    if mode == "fast":
        samples = CNVpos

    q = 0
    for name in samples:
        abAmp = aberrantAmpliconsPerSample(name, reads_norm, CNVneg, conta=conta)
        if abAmp.shape != (0, ):
            genes = np.unique([i.split('_')[0] for i in abAmp])
            for gene in genes:
                r = re.compile(gene)
                abEx = list(filter(r.match, abAmp))
                tmp = reads.filter(regex=f"^{gene}", axis=0)

                score = scoreAmplif(len(abEx), tmp.shape[0], reads.shape[0])
                amplif = amplifEvalGene(reads_norm, gene, CNVneg, name)

                if score > scoreThreshold:
                    f.loc[q] = [run, name, gene, amplif, score]
                    q = q + 1
    if verbose:
        print(f"{str(f.shape[0])} aberrant regions found in {str(len(np.unique(f['Sample name'])))} samples.\n")

    return f


###########################################
#               MAIN                      #
###########################################
def main(args):
    if args.lib_ressources is None:
        lib_ressources = f"{os.path.dirname(os.path.abspath(__file__))}/ressources"
    else:
        lib_ressources = args.lib_ressources
    # Create or open the reads matrix
    if args.skip is None:
        reads = createReadsMatrix(pathToBam=args.input, bedFile=args.bed, verbose=args.verbose, output=args.readsMatrixOuptut)
    else:
        reads = pd.read_csv(args.skip, sep="\t", index_col=0)

    # Filter the reads matrix
    filteredReads, filteredS, filteredT = filterReads(reads=reads, N=args.minReads, regtar=args.regTargets, regsamp=args.regSample)

    if args.verbose:
        if len(filteredS) > 0:
            print("Filtered samples:")
            for i in filteredS:
                print(i)

    # Normalize the reads matrix
    normReads = normalizeReads(filteredReads)

    # Find the aberrant samples
    CNVpos, CNVneg = aberrantSamples(filteredReads, conta=args.contaSamples, verbose=args.verbose)

    # Find the aberrant targets
    final = aberrantAmpliconsFinal(filteredReads, normReads, CNVpos, CNVneg, scoreThreshold=args.scoreThreshold, conta=args.contaTargets, mode=args.mode, run=args.run, verbose=args.verbose)

    if args.save:
        if not os.path.isdir(args.output):
            cmd = ["mkdir", args.output]
            subprocess.check_output(cmd)
        final.to_csv(f"{args.output}/{args.run}.tsv", sep="\t")

    # Generate the report
    generateReport(final, output_dir=args.output, reads=normReads, ressources=lib_ressources, CNVneg=CNVneg)

    if args.verbose:
        print("ifCNV analysis done succesfully !\n")

    if args.autoOpen:
        try:
            os.startfile(f"{args.output}/run.html")
        except AttributeError:
            try:
                subprocess.call(['open', f"{args.output}/run.html"])
            except:
                print('Could not open URL')
