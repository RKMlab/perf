/*-- BAR PLOT -------------------------------------------------------------------------------*/
const barDatum = function(data, n, s, repeats, rs) {
    let datum = []
    let orderedRepeats = _.flatMap(repeatSet, function(d) { return d });
    if (repeats) {
        for (let d in orderedRepeats) { d = orderedRepeats[d]; if (repeats.indexOf(d) > -1) { datum.push([d, _.sum(data[d])]); } }
        if (rs == 1) {
            datum = _.sortBy(datum, function(d) { return d[1]; });
            datum = datum.reverse();
        }
    } else if (n) {
        for (let d in data) { datum.push([d, _.sum(data[d])]) };
        datum = _.sortBy(datum, function(d) { return d[1]; });
        if (!(s)) { datum = datum.reverse(); }
        datum = datum.slice(0, n);
    }
    return datum;
}

const barPlot = function(selectionType, repeats) {
    let plotData;
    if (selectionType == 'sort') {
        const numBars = document.getElementById('bar-num').value;
        const sort = document.getElementById('bar-select').selectedIndex;
        plotData = barDatum(plotInfo['len'], numBars, sort);
    } else if (selectionType == 'repeat') {
        const barSortRadios = document.getElementsByName('bar-sort');
        let repeatSort;
        for (let i in barSortRadios) {
            let radio = barSortRadios[i];
            if (radio.tagName == 'INPUT') { if (radio.checked) { repeatSort = radio.value; } }
        }
        plotData = barDatum(plotInfo['len'], 0, 0, repeats, repeatSort);
    }

    Highcharts.chart('bar-plot-svg', {
        chart: { type: 'column', marginTop: 40, marginLeft: 100 },
        title: { text: null },
        xAxis: { type: 'category', title: { text: 'Repeat Class' }, labels: { rotation: -45, style: { fontSize: '12px', fontFamily: 'Verdana, sans-serif' } } },
        yAxis: { min: 0, title: { text: 'Frequency' } },
        legend: { enabled: false },
        tooltip: { pointFormat: 'Frequency: <b>{point.y}</b>' },
        series: [{ name: 'Frequency', data: plotData, dataLabels: { enabled: false, } }]
    });
}

document.getElementById('bar-num').oninput = function() { barPlot('sort') };
document.getElementById('bar-select').onchange = function() { barPlot('sort') };
let barSetValues = [];
let barPlotRepeats;
let barRepeatSelect = new SlimSelect({ select: "#bar-repeats-sel", placeholder: 'Select Repeats', });

barRepeatSelect.beforeOnChange = function(info) {
    if (info.length != 0) {
        let allValues = _.map(info, 'value');
        barSetValues = allValues;
        let lastValue = info[info.length - 1].value;
        if (lastValue == 'select-all') {
            for (let k in repeatSet) {
                barSetValues = barSetValues.concat([`select-all-${k}`])
                barSetValues = barSetValues.concat(repeatSet[k]);
                barSetValues = _.uniq(barSetValues);
                barRepeatSelect.set([]);
                barRepeatSelect.set(barSetValues);
            }
        } else if (lastValue.length > 10 && lastValue.slice(0, 10) == 'select-all') {
            let kmer = lastValue.slice(11, lastValue.length);
            barSetValues = allValues.concat(repeatSet[kmer]);
            barSetValues = _.uniq(barSetValues);
            barRepeatSelect.set([]);
            barRepeatSelect.set(barSetValues);
        }
        barPlotRepeats = _.filter(barSetValues, function(d) {
            return d.slice(0, 10) != 'select-all';
        });
        barPlot('repeat', barPlotRepeats);
    }
}

barRepeatSelect.onChange = function(info) {
    let currentValues = _.map(info, 'value');
    if (currentValues.length == 0 && barSetValues.length - currentValues.length == 1 && barSetValues[0].slice(0, 10) == 'select-all') {
        // pass;
    } else if (barSetValues.length - currentValues.length == 1) {
        let removedValue = (_.difference(barSetValues, currentValues))[0];
        if (removedValue == 'select-all') {
            currentValues = [];
            barRepeatSelect.set([]);
            barRepeatSelect.set(currentValues);
            barSetValues = currentValues;
            barPlotRepeats = _.filter(barSetValues, function(d) {
                return d.slice(0, 10) != 'select-all';
            });
            barPlot('repeat', barPlotRepeats);
        } else if (removedValue.length > 10 && removedValue.slice(0, 10) == 'select-all') {
            let kmer = removedValue.slice(11, removedValue.length);
            let tempCurrentValues = _.difference(currentValues, repeatSet[kmer]);
            currentValues = tempCurrentValues;
            barSetValues = currentValues;
            if (_.intersection(repeatSet[kmer], barPlotRepeats).length > 0) {
                barRepeatSelect.set([]);
                barRepeatSelect.set(currentValues);
                barPlotRepeats = _.filter(barSetValues, function(d) {
                    return d.slice(0, 10) != 'select-all';
                });
                barPlot('repeat', barPlotRepeats);
            }
        } else {
            barSetValues = currentValues;
            barPlotRepeats = _.filter(barSetValues, function(d) {
                return d.slice(0, 10) != 'select-all';
            });
            barPlot('repeat', barPlotRepeats);
        }
    } else if (currentValues.length == 0) {
        //pass;
    }
}

const barSortRadios = document.getElementsByName('bar-sort');
for (let i in barSortRadios) {
    let radio = barSortRadios[i];
    if (radio.tagName == 'INPUT') {
        radio.onchange = function() { barPlot('repeat', barPlotRepeats); }
    }
}


let initial = 1;
const placeholderBarPlot = function() {
    if (initial) {
        barPlot('sort');
        initial = 0
    }
}