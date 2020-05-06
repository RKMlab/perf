/*
    The "main.js" contains core JS component supporting the PERF-Analysis module.

    This web application is developed with Semantic-ui frame work.
    The charts are build using Apex-Chart js charting library.

    All the data for the report is derived from analyse_data.js
    data = {info: {seqInfo: {}, repInfo: {}, plotInfo: {}}}

    plotInfo i s a dictionary with key as the repeat class and value as a dictionary
    plotInfo: { REPEAT_CLASS: { LENGTH: FREQUENCY } }

*/

// Updating report data
for (const key in data.info.seqInfo){$(`.value.${key}`).html(data.info.seqInfo[key])};
for (const key in data.info.repInfo){if (key != 'lenFrequency') { $(`.value.${key}`).html(data.info.repInfo[key]); }};

const menuLayout = function(){
    const w = window.innerWidth;
    if (w < 800) { 
        const navmenu = document.getElementById('navmenu');
        navmenu.classList.remove('vertical');
        navmenu.parentElement.style.width = '100%';
        document.getElementById('content-display').style.width = '100%';
    }
    
    else {
        const navmenu = document.getElementById('navmenu')
        navmenu.classList.add('vertical');
        navmenu.parentElement.style.width = '5%';
        document.getElementById('content-display').style.width = '95%';
    }
}
window.onresize = function(){ menuLayout(); }

const numPrefixObj = ["Monomer","Dimer","Trimer","Tetramer","Pentamer","Hexamer","Heptamer","Octamer","Nonamer","Decamer","Undecamer","Dodecamer","Tridecamer","Tetradecamer","Pentadecamer","Hexadecamer","Heptadecamer","Octadecamer","Nonadecamer","Icosamer","Uncosamer","Docosamer","Tricosamer","Tetracosamer","Pentacosamer","Hexacosamer","Heptacosamer","Octacosamer","Nonacosamer","Triacontamer","Untriacontamer","Dotriacontamer","Tritriacontamer","Tetratriacontamer","Pentatriacontamer","Hexatriacontamer","Heptatriacontamer","Octatriacontamer","Nonatriacontamer","Tetracontamer","Untetracontamer","Dotetracontamer","Tritetracontamer","Tetratetracontamer","Pentatetracontamer","Hexatetracontamer","Heptatetracontamer","Octatetracontamer","Nonatetracontamer","Pentacontamer"]
const plotData = data.info.repInfo.lenFrequency;
const allRepClasses = data.info.repInfo.allRepClasses;

$('.ui.dropdown').dropdown();
$('.chart .item').tab();
$('.anno-chart .item').tab();
$('.ui .units').dropdown({values: [{name: 'length', value: 1, selected:true}, {name: 'units', value: 0}]});

/* 
    Bar graph
    - For bar graph we curate the data in barData with Repeat class as the key and [bases, frequency] as the value.

    - bar_activeSelected is the list of Repeat classes which are considered to be plotted.
        - can be selected by sort selection which is bar_sortSelected.
        - or by repeat selection dropdown stored in bar_repSelected.
        - sorted_barKeys stores the all repClasses sorted based on the datatype selected.
    
    The dataflow for the barChart is as follows.
    - There are two plot buttons dedicated individually to either plot with data based sorted keys or desired set of keys.
    - Based on which plot button is pressed that repeat set gets placed in the bar_activeSelected variable.
    - And then subsequently the data for the bar_activeSelected keys is plotted.
    - The data type and the sort customisation of the plots only deal with the bar_activeSelected keys.

*/

const barData = {};
allRepClasses.forEach(function(e){
    if ((Object.keys(plotData).indexOf(e) != -1) && (plotData[e] != 0)) {
        const lengths = _.map(Object.keys(plotData[e]), d => { return parseInt(d); });
        let frequency = 0; let bases = 0;
        for (let l in lengths) { 
            l = lengths[l]; frequency += _.sum(plotData[e][l]); bases += _.sum(plotData[e][l]) * l;
        }
        barData[e] = [bases, frequency];
    } else { barData[e] = [0, 0]; }
})
let bar_dataType = 1;
let bar_sortSelected = [];
let bar_sortOrder = 1;
let bar_numReps = 10;
let bar_repSelected = [];
let sorted_barKeys = _.sortBy(Object.keys(barData), k => { return barData[k][bar_dataType]; });
let bar_activeSelected = [];

const bar_options = {
    chart: { type: 'bar' },
    plotOptions: { bar: { horizontal: false, columnwidth: '55%' } },
    series: [{ data: [] }],
    dataLabels: { enabled: false },
    yaxis: { 'title': { 'text': 'Frequency', 'style': { 'fontSize': '16px', 'font-weight': 'bold' } }},
    xaxis: { categories: [], 'title': { 'text': 'Repeat Class', 'style': { 'fontSize': '16px', 'font-weight': 'bold' } }},
    title: { text: 'Repeat Frequency', align: 'left' }
}
const bar_chart = new ApexCharts(document.querySelector('#bar-plot-area'), bar_options);
bar_chart.render();
const plotBar = function(keys){ 
    const values = [];
    keys.forEach(function(e){ values.push(barData[e][bar_dataType]); });
    const name = (bar_dataType == 1) ? 'Frequency' : 'Bases';
    bar_chart.updateOptions({series: [{'name': name, data: values}], yaxis: { title: { 'text': name, 'style': { 'fontSize': '16px', 'font-weight': 'bold' } } }, xaxis: {categories: keys}, animate: true})
}

$('#bar-numRep').change(function(){ bar_numReps = this.value; });
$('.ui .dropdown.sort-order').dropdown({
    values: [{name: 'top', value: 1, selected:true}, {name: 'bottom', value: 0}],
    onChange: function(value) { bar_sortOrder = value; }
});
$('#bar-sortPlot-button').click(function(){ 
    if (bar_sortOrder == 1) { bar_sortSelected = sorted_barKeys.slice(sorted_barKeys.length - bar_numReps); bar_sortSelected.reverse();}
    else { bar_sortSelected = sorted_barKeys.slice(0, bar_numReps); }
    bar_activeSelected = bar_sortSelected; plotBar(bar_activeSelected); 
})
$('#bar-sortPlot-button').trigger("click")

$("#bar-repeat-select").multiSelect({
    selectableOptgroup: true,
    afterSelect: function(d){ d.forEach(function(e){ if (bar_repSelected.indexOf(e) == -1) { bar_repSelected.push(e) } })},
    afterDeselect: function(d){ d.forEach(element => { bar_repSelected.splice(bar_repSelected.indexOf(element), 1); }); } 
});
$('#bar-repPlot-button').click(function(){ 
    bar_repSelected = _.sortBy(bar_repSelected, o => {return allRepClasses.indexOf(o)}); 
    bar_activeSelected = bar_repSelected; plotBar(bar_activeSelected); 
})

// Once the data type is selected the global variable bar_dataType changes
// also the sorted_barKeys is updated based on the datatype.
$('.ui.checkbox.bar').checkbox({ onChange: function(val){
    bar_dataType = this.value;
    sorted_barKeys = _.sortBy(Object.keys(barData), k => { return barData[k][bar_dataType]; }); 
    plotBar(bar_activeSelected)
}});

$('#asort-alpha').click(function(){ bar_activeSelected = bar_activeSelected.sort(function(a, b){ return allRepClasses.indexOf(a) - allRepClasses.indexOf(b) }); plotBar(bar_activeSelected); })
$('#dsort-alpha').click(function(){ bar_activeSelected = bar_activeSelected.sort(function(a, b){ return allRepClasses.indexOf(a) - allRepClasses.indexOf(b) }); bar_activeSelected.reverse(); plotBar(bar_activeSelected); })
$('#asort-num').click(function(){ bar_activeSelected = _.sortBy( bar_activeSelected, k => { return barData[k][bar_dataType];}); plotBar(bar_activeSelected); })
$('#dsort-num').click(function(){ bar_activeSelected = _.sortBy( bar_activeSelected, k => { return barData[k][bar_dataType];}); bar_activeSelected.reverse(); plotBar(bar_activeSelected); })


/* 
    Pie graph
    - For pie chart as we deal with frequency and bases data of repeat classes we continue using barData for data retrieval.

    - pie_activeSelected is the list of Repeat classes which are considered to be plotted.
        - As there is only way to select repeats there will be no updating the repeats.
    
    The dataflow for the barChart is as follows.
    - There are two plot buttons dedicated individually to either plot with data based sorted keys or desired set of keys.
    - Based on which plot button is pressed that repeat set gets placed in the bar_activeSelected variable.
    - And then subsequently the data for the bar_activeSelected keys is plotted.
    - The data type and the sort customisation of the plots only deal with the bar_activeSelected keys.

*/

let pie_activeSelected = allRepClasses;
let pie_dataType = 1;
let pie_group = true;
$("#pie-kmer-toggle").checkbox({
    onChecked: function(){ pie_group = true; plotPie(pie_activeSelected) },
    onUnchecked: function(){ pie_group= false; plotPie(pie_activeSelected) }
});
$("#pie-repeat-select").multiSelect({
    selectableOptgroup: true,
    afterSelect: function(d){ d.forEach(function(e){ if (pie_activeSelected.indexOf(e) == -1) { pie_activeSelected.push(e) } })},
    afterDeselect: function(d){ d.forEach(element => { pie_activeSelected.splice(pie_activeSelected.indexOf(element), 1); }); } 
});
$(".ui.checkbox.pie.radio.pie-data-type").checkbox({ onChange: function(){ pie_dataType = this.value; plotPie(pie_activeSelected) }});
const pie_options = {
    chart: { type: 'pie' },
    labels: ['Monomers', 'Dimers', 'Trimers', 'Tetramers', 'Pentamers', 'Hexamers'],
    series: [10, 10, 10, 10, 10, 10],
    responsive: [{ breakpoint: 480 }],
    colors: ["#3366cc", "#dc3912", "#ff9900", "#109618", "#990099", "#0099c6", "#dd4477", "#66aa00", "#b82e2e", "#316395", "#994499", "#22aa99", "#aaaa11", "#6633cc", "#e67300", "#8b0707", "#651067", "#329262", "#5574a6", "#3b3eac"],
    // theme: { monochrome: { enabled: true }}
}
let pie_chart = new ApexCharts( document.querySelector("#pie-plot-area"), pie_options );
pie_chart.render();

const plotPie = function(keys){ 
    let values = [];
    keys = keys.sort(function(a, b){ return allRepClasses.indexOf(a) - allRepClasses.indexOf(b) })
    keys.forEach(function(e){ values.push(barData[e][pie_dataType]); });
    if (pie_group == true) { 
        values = [];
        let group_keys = [];
        const kmer_lengths = _.uniq(_.map(keys, e => { return e.length; })).sort();
        kmer_lengths.forEach( e => { group_keys.push(numPrefixObj[e - 1]); values.push(0); });
        keys.forEach(e => { values[kmer_lengths.indexOf(e.length)] += barData[e][pie_dataType]; })
        keys = group_keys;
    }
    pie_chart.updateOptions({labels:keys, series: values, animate: true})
}

$("#pie-plot-button").click(function(){
    pie_activeSelected = _.sortBy(pie_activeSelected, o => {return allRepClasses.indexOf(o)});
    plotPie(pie_activeSelected);
});

//  Initialising Pie Chart
plotPie(pie_activeSelected);

/* 
    Line graph
    - We retrieve the data for line plots from plotData.

    - Line plot has options to
        - Select repeats which will be saved in line_activeSelected variable
        - Length range in which it has to plot the data
    
    - Required data
        - The minimum length/units criteria is retrieved from the repInfo.

    The dataflow for the lineChart is as follows.
    - The data flow is prestty simple from the options selected the relavant data is retrieved and plotted

*/

let minLength = data.info.repInfo.minLength;
let minUnits = data.info.repInfo.minUnits;
let minRange = 12;
let maxRange = 50;
let line_dataType = 1;
let line_activeSelected = ['A', 'C'];


$('.ui .dropdown.units').dropdown({
    values: [{name: 'length', value: 1, selected:true}, {name: 'units', value: 0}],
    onChange: function(value) { line_dataType = value;}
});
$("#line-repeat-select").multiSelect({
    selectableOptgroup: true,
    afterSelect: function(d){ d.forEach(function(e){ if (line_activeSelected.indexOf(e) == -1) { line_activeSelected.push(e) } }) },
    afterDeselect: function(d){ d.forEach(element => { line_activeSelected.splice(line_activeSelected.indexOf(element), 1); }); } 
});
$('#line-min-len').change(function(){ minRange = parseInt(this.value); plotLine(line_activeSelected); });
$('#line-max-len').change(function(){ maxRange = parseInt(this.value); plotLine(line_activeSelected); });
const line_options = {
    chart: { type: 'line', zoom: { enabled: false }},
    dataLabels: { enabled: false },
    stroke: { curve: 'straight', width: 2 },
    series: [],
    title: { text: 'Repeat sequence length(bp) vs Abundance', align: 'left' },
    grid: { row: { colors: ['#f3f3f3', 'transparent'], opacity: 0.5 } },
    tooltip: { x: {
        formatter: function(val) { return `Length: ${val}bp` }
    }},
    markers: {size: 0},
    yaxis: { title: { text: 'Frequency', 
             style: { 'fontSize': '16px', 'font-weight': 'bold' } }},
    xaxis: { title: { text: 'Length (bp)', 
             style: { 'fontSize': '16px', 'font-weight': 'bold' } },
            //  labels: { format: '%d' },
             tickAmount: parseInt((maxRange - minRange)/2) },
    legend: { position: 'top' }
}
const line_chart = new ApexCharts( document.querySelector("#line-plot-area"), line_options );
line_chart.render();

const plotLine = function(keys) {
    const xValues = _.range(minRange, maxRange + 1);
    const series = [];
    keys.forEach(function(key){
        const data = [];
        if (line_dataType == 0) {
            for (let i = minRange; i <= maxRange; i++) {
                let val = 0;
                for ( let j = 0; j < key.length; j++) { const repLen = (i*key.length) + j; const v = (plotData[key][repLen]) ? _.sum(plotData[key][repLen]) : 0; val += v; }
                data.push(val);
            }
        }
        else { for (let i = minRange; i <= maxRange; i++){ const val = (plotData[key][i]) ? _.sum(plotData[key][i]) : 0; data.push(val); } }
        series.push({'name': key, 'data': data});
    })
    line_chart.updateOptions({series: series, xaxis: {categories: xValues}})
}
$('#line-plot-button').click(function(){ plotLine(line_activeSelected); })

// Initialsing line chart
plotLine(line_activeSelected);
