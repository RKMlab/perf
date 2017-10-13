const sizeSuffix = function(size) {
    if (size < 1000) { return size + 'bp'; } else if (size < 1000000) { return (size / 1000).toFixed(2) + 'Kb'; } else if (size < 1000000000) { return (size / 1000000).toFixed(2) + 'Mb'; } else { return (size / 1000000000).toFixed(2) + 'Gb'; }
}

let numPrefixObj = { 1: 'Monomer', 2: 'Dimer', 3: 'Trimer', 4: 'Tetramer', 5: 'Pentamer', 6: 'Hexamer', 7: 'Septamer', 8: 'Octamer', 9: 'Enneamer', 10: 'Decamer' }
let plotInfo = data['info']['plotInfo'];
let minLength = parseInt(data['info']['minLength']);
let minUnits = parseInt(data['info']['minUnits']);
let kmerLengths = _.uniq(_.map(Object.keys(plotInfo['len']), function(o) { return o.length; })).sort();
let kMers = _.map(kmerLengths, function(o) { if (_.map(Object.keys(numPrefixObj), parseFloat).indexOf(o) != -1) { return numPrefixObj[o] } else { return o + '-kmer' } });
let kmerObj = {};
for (let k in kMers) { kmerObj[kmerLengths[k]] = kMers[k]; }
let repeatSet = {};
for (let k in kMers) { repeatSet[kMers[k]] = []; }
for (let repeat in plotInfo['len']) { repeatSet[kmerObj[repeat.length]].push(repeat); }
for (let group in repeatSet) { repeatSet[group].sort(); }
let basicInfoKeys = ['name', 'genomeSize', 'GC', 'numSeq', 'numRepClass', 'totalRepBases', 'totalRepFreq', 'percentGenomeCovered', 'repDensityByFreq', 'repDensityByBases'];
for (let k in basicInfoKeys) {
    k = basicInfoKeys[k];
    let spanElement = document.getElementById(k);
    if (k == 'totalRepBases' || k == 'genomeSize') { spanElement.innerHTML = sizeSuffix(data['info'][k]); } else { spanElement.innerHTML = data['info'][k]; }
    spanElement.setAttribute('title', data['info'][k]);
}

const populateTable = function(id, data) {
    const tableBodyElement = document.querySelector("#" + id + "> tbody");
    const keys = ['seq', 'start', 'end', 'repClass', 'repLength', 'repOri', 'repUnit', 'actualRep']
    const keyObj = { seq: "Sequence id", start: "Start", end: "End", repClass: "Repeat class", repLength: "Length", repOri: "Strand", repUnit: "Repeat units", actualRep: "Actual repeat" }
    for (let o in data) {
        let obj = data[o];
        let rowElement = document.createElement("tr");
        for (let key in keys) {
            key = keys[key];
            let cellElement = document.createElement("td");
            cellElement.setAttribute('data-label', keyObj[key]);
            cellElement.setAttribute('title', obj[key]);
            cellElement.innerHTML = obj[key];
            rowElement.appendChild(cellElement)
        }
        tableBodyElement.appendChild(rowElement);
    }
}

populateTable('most-repeat-units-table', data['info']['mostRepeatUnits']);
populateTable('longest-repeats-table', data['info']['longestRepeats']);
const addOptionsToMultiSelect = function(id) {
    let multiSelect = document.getElementById(id);
    let selectAll = document.createElement("option");
    selectAll.setAttribute("value", "select-all");
    selectAll.setAttribute("id", "select-all");
    selectAll.innerHTML = "Select All";
    multiSelect.appendChild(selectAll);
    for (let kmer in kMers) {
        kmer = kMers[kmer];
        let kmerOptGroup = document.createElement("optgroup");
        kmerOptGroup.setAttribute("label", kmer);
        let selectAll = document.createElement("option");
        selectAll.setAttribute("value", `select-all-${kmer}`)
        selectAll.innerHTML = 'Select All ' + kmer;
        kmerOptGroup.appendChild(selectAll);
        for (let repeat in repeatSet[kmer]) {
            repeat = repeatSet[kmer][repeat];
            let repeatOpt = document.createElement("option");
            repeatOpt.setAttribute("value", repeat);
            repeatOpt.innerHTML = repeat;
            kmerOptGroup.appendChild(repeatOpt);
        }
        multiSelect.appendChild(kmerOptGroup);
    }
}

addOptionsToMultiSelect('bar-repeats-sel');
addOptionsToMultiSelect('pie-repeats-sel');
addOptionsToMultiSelect('line-repeats-sel');