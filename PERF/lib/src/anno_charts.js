const repeatAnnoDist = function(data, repeats, percent=true) {
    const annoKeyObj = {
        'Exon': ['EP', 'EN'],
        'Intron': ['IP', 'IN'],
        'Genic': ['GP', 'GN'],
        'Intergenic': ['DP', 'DN'],
    };
    const promKeyObj = {
        'Promoter': ['EP', 'IP', 'GP', 'DP'],
        'Non-Promoter': ['EN', 'IN', 'GN', 'DN']
    };

    let annoKeys = Object.keys(annoKeyObj);
    let promKeys = Object.keys(promKeyObj);
 
    let annodata = {};
    if (repeats === "all") { 
        let obj = {};
        let annototal = 0;
        let promtotal = 0;
        annoKeys.forEach(a => {
            let val = 0;
            Object.keys(data).forEach( r => {
                annoKeyObj[a].forEach( a => {
                    val += data[r][a]
                })
            })
            obj[a] = val;
            annototal += val;
        });
        
        promKeys.forEach(a => {
            let val = 0;
            Object.keys(data).forEach( (r, i) => {
                promKeyObj[a].forEach( a => { val += data[r][a] })
            })
            obj[a] = val;
            promtotal += val;
        });
        
        if (percent){
            annoKeys.forEach(a => { obj[a] = ((obj[a]*100)/annototal).toFixed(2) })
            promKeys.forEach(a => { obj[a] = ((obj[a]*100)/promtotal).toFixed(2) })
        }
        annodata = obj
    }
    
    else {
        
        for (let rep of repeats) {
            const repindex = Object.keys(data).indexOf(rep);
            if (repindex > -1) {
              let obj = {};
              annoKeys.forEach(a => {
                  obj[a] = _.sum(annoKeyObj[a].map(d => { return data[rep][d]; }));
              });
              promKeys.forEach(p => {
                  obj[p] = _.sum(promKeyObj[p].map(d => { return data[rep][d]; }));
              });
              annodata[rep] = obj;
            }
        }
        for (let rep of repeats) {
            const repindex = Object.keys(data).indexOf(rep);
            if (repindex > -1) {
              const annototal = _.sum(annoKeys.map(a => {return annodata[rep][a] }));
              const promtotal = _.sum(promKeys.map(p => {return annodata[rep][p] }));
              if (percent) {
                  for (let a of annoKeys) { annodata[rep][a] = ((annodata[rep][a]*100)/annototal).toFixed(2); }
                  for (let a of promKeys) { annodata[rep][a] = ((annodata[rep][a]*100)/promtotal).toFixed(2); }
              }
            }
        }
    }


    return annodata;
}


const kmerAnnoDist = function(plotdata, freqdata, stacktype="annotation"){
    const annoKeys = ['Exon', 'Intron', 'Intergenic', 'Genic', 'Promoter', 'Non-Promoter'];
    const repeatsObj = {}
    const kmerFreq = {}
    const kmers = [];
    repeats.forEach(d => {
      repeatsObj[d.kmer] = _.map(d.repeats, 'class'); 
      kmerFreq[d.kmer]
      kmers.push(d.kmer);
    })
    const annoKmerData = {};
    for (const anno of annoKeys) {
      annoKmerData[anno] = {'Monomer': 0, 'Dimer': 0, 'Trimer': 0, 'Tetramer': 0, 'Pentamer': 0, 'Hexamer': 0};
    }
    
    for (let kmer in repeatsObj) {
      const classes = repeatsObj[kmer];
      const kmerData = repeatAnnoDist(plotdata, classes, false)
      const repFreqData = repeatFrequency(freqdata, classes, 'kmer', 'freq');
      kmerFreq[kmer] = _.sumBy(repFreqData, 'value');
      for (let rep in kmerData) { 
        for (let anno in kmerData[rep]){
          annoKmerData[anno][kmer] += kmerData[rep][anno]
        } 
      }
    }

    const outdata = {data: {} };
    if (stacktype === 'kmer') {
      outdata.keys = kmers;
      for (const anno in annoKmerData) {
        outdata.data[anno] = [];
        for (const kmer of kmers) {
          // const total = kmerFreq[kmer];
          const annototal = annoKmerData['Exon'][kmer] + annoKmerData['Intron'][kmer] + annoKmerData['Intergenic'][kmer] + annoKmerData['Genic'][kmer]
          if (annototal > 0) {
            outdata.data[anno].push(((annoKmerData[anno][kmer]*100)/annototal).toFixed(2))
          }
        }
      }
    }
    else {
      for (const kmer of kmers) {
          outdata.data[kmer] = []
          outdata.keys = []
          for (const anno in annoKmerData) {
            const total = _.sum(Object.values(annoKmerData[anno]));
            if (total > 0) {
              outdata.keys.push(anno);
              outdata.data[kmer].push(((annoKmerData[anno][kmer]*100)/total).toFixed(2));
            }
          }
      }
    }
    return outdata;
}

const tssHistData = function(data, bins, repeats, bin, datatype) {
    const values = _.range(-4975, 4976, 50);

    const stepSize = parseInt(bin/50);
    const steps = parseInt((values.length)/stepSize);

    const start = -5000 + parseInt(bin/2);
    const end = 5001 - parseInt(bin/2);
    const binCenters = _.range(start, end, bin);

    const y = {};
    repeats.forEach(r => {
        y[r] = Array(binCenters.length).fill(0);
        const repindex = Object.keys(data).indexOf(r);
        if (repindex > -1) {
            for (let i=0; i<steps; i++){
                let val = _.sum(data[r].slice(i*stepSize, (i+1)*stepSize));
                y[r][i] = val;
            }
        }

    })

    for (let rep of repeats) {
        if (datatype === 'density'){
            const total = _.sum(y[rep]);
            let a;
            if (total != 0) {
                a = _.map(y[rep], (d, i) => { 
                    let val = (d)/total; 
                    return [binCenters[i], val.toFixed(3)]; 
                });
            }
            else {
                a = _.map(y[rep], (d, i) => { 
                    return [binCenters[i], 0]; 
                });
            }
            y[rep] = a;
        }
        else {
            let a = _.map(y[rep], (d, i) => {
                 let val = d; return [binCenters[i], val]; 
            })
            y[rep] = a;
        }
        y[rep].push([5000, y[rep][y[rep].length - 1][1]]);
        y[rep].unshift([-5000, y[rep][0][1]]);
    }

    return { data: y };
}

if (data.info.annoInfo) { 
    
    const annostackbar_activeSelected = ['A', 'C']; //allRepClasses;
    let stack_group = false;
    
    $("#anno-stackbar-repeat-select").multiSelect({
        selectableOptgroup: true,
        afterSelect: function(d){ d.forEach(function(e){ if (annostackbar_activeSelected.indexOf(e) == -1) { annostackbar_activeSelected.push(e) } })},
        afterDeselect: function(d){ d.forEach(element => { annostackbar_activeSelected.splice(annostackbar_activeSelected.indexOf(element), 1); }); } 
    });
    
    const annostackbar_data = repeatAnnoDist(data.info.annoInfo.repAnno, allRepClasses, false);
    var annostackbar_options = {
        series: [],
        chart: { type: 'bar', stacked: true, stackType: '100%' },
        plotOptions: {},
        stroke: { width: 1, colors: ['#fff'] },
        title: { text: 'Repeat Genomic distribution' },
        xaxis: {  },
        tooltip: { y: { formatter: function (val) { return val } } },
        fill: { opacity: 1 },
        legend: { position: 'top', horizontalAlign: 'left', offsetX: 40 }
    };
    var annostackbar_chart = new ApexCharts(document.querySelector("#anno-stackbar-plot-area"), annostackbar_options);
    annostackbar_chart.render();
    
    const annoLabels = ['Exon', 'Intron', 'Genic', 'Intergenic']
    const plot_annostackbar = function(){
        const annostackbar_series = []
        if (!(stack_group)) {
            annoLabels.forEach(function(a){ 
                const d = _.map(annostackbar_activeSelected, o => {
                    return parseFloat(annostackbar_data[o][a])
                })
                annostackbar_series.push({name: a, data: d});
            })
            annostackbar_chart.updateOptions({series: annostackbar_series, xaxis: {categories: annostackbar_activeSelected}});
        }

        else {
            annoLabels.forEach(function(a){ 
                const d = _.map(annostackbar_activeSelected, o => {
                    return parseFloat(annostackbar_data[o][a])
                })
                annostackbar_series.push({name: a, data: [_.sum(d)]});
            })
            annostackbar_chart.updateOptions({series: annostackbar_series, xaxis: {categories: ['All selected repeats']}});
        }
    }

    $('.ui.checkbox.anno-stackbar').checkbox({ onChange: function(){
        stack_group = !(stack_group);
        plot_annostackbar();
    }});

    $("#anno-stackbar-plot-button").click(function(){ plot_annostackbar(); });
    plot_annostackbar();


    const annoarea_activeSelected = ['A', 'C']; //allRepClasses;
    let binSize = 500;
    $("#anno-area-repeat-select").multiSelect({
        selectableOptgroup: true,
        afterSelect: function(d){ d.forEach(function(e){ if (annoarea_activeSelected.indexOf(e) == -1) { annoarea_activeSelected.push(e) } })},
        afterDeselect: function(d){ d.forEach(element => { annoarea_activeSelected.splice(annoarea_activeSelected.indexOf(element), 1); }); } 
    });

    $('.ui .dropdown.bin-size').dropdown({
        values: [
            {name: 100, value: 100},
            {name: 200, value: 200},
            {name: 500, value: 500, selected:true},
            {name: 1000, value: 1000}
        ],
        onChange: function(value) { binSize = value; }
    });

    var annoarea_options = {
        series: [],
        chart: { type: 'area'},
        plotOptions: {},
        stroke: { width: 1 },
        title: { text: 'Average repeat distribution around TSS' },
        xaxis: { type: 'numeric', min: -5000, max: 5000, tickAmount: 10000/binSize, axisTicks: { height: 8 }},
        tooltip: { 
            y: { formatter: function (val) { return val } },
            x: { formatter: function (val) { return `${val-parseInt(binSize/2)}bp - ${val+parseInt(binSize/2)}bp` } } 
        },
        fill: { opacity: 1 },
        legend: { position: 'top', horizontalAlign: 'left', offsetX: 40 }
    };
    var annoarea_chart = new ApexCharts(document.querySelector("#anno-area-plot-area"), annoarea_options);
    annoarea_chart.render();

    const plot_annoarea = function() {
        const annoarea_data = tssHistData(data.info.annoInfo.TSS_dist, data.info.annoInfo.TSS_histBinEdges, annoarea_activeSelected, binSize, 'density')['data'];
        const series = []
        for (let rep of Object.keys(annoarea_data)) { series.push({ name: rep, data: _.map(annoarea_data[rep], o => { return [o[0], parseFloat(o[1])];})}) }
                annoarea_chart.updateOptions({ series: series, xaxis: { type: 'numeric', min: -5000, max: 5000, tickAmount: 20, axisTicks: { height: 8 }} });
    }

    $("#anno-area-plot-button").click(function(){ plot_annoarea(); });
    plot_annoarea();

}

else { document.getElementById('anno-charts-main').style.display = 'none'; }