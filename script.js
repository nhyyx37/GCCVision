// Main configuration parameters
let config = {
    margin: { top: 60, right: 200, bottom: 60, left: 100 }, // Increase right margin to accommodate proportion info
    chromosomeHeight: 25,
    chromosomeSpacing: 50,     
    maxWidth: 1000,
    siteRadius: 2,
    colors: {
        'Hom_A': '#e74c3c',
        'Hom_B': '#3498db',
        'Het_AB': '#2ecc71'
    },
    histogramBinSize: 100000,  
    geneConnectorOffset: 2,    
    displayMode: 'histogram',  // Default to histogram mode
    backgroundLabelOffset: 20,  // Offset for parental background proportion labels
    geneFontSize: 10,  // Gene name font size, default 10px
    backgroundFontSize: 10,  // Genetic background proportion font size, default 10px
    chrBorderWidth: 1,  // Add default chromosome border width
    geneColorMode: 'auto',  // Gene color display mode: 'custom'(use specified color), 'auto'(calculate based on sites), 'black'(all black)
    geneAngle: 0,  // Gene name display angle: 0 degrees for horizontal (all above), 45 degrees for 45-degree angle (alternating above/below)
    alternateSides: false,  // Control whether gene names alternate above/below, default no alternation
    labels: {
        'Hom_A': 'Hom_A',
        'Hom_B': 'Hom_B',
        'Het_AB': 'Het_AB',
        'ParentA': 'Parent A:',
        'ParentB': 'Parent B:'
    }
};

// Default data
const defaultChrData = 
`Chr01	56831624
Chr02	48577505`;

const defaultSiteData = 
`Chr01	1959	Hom_A
Chr01	2528	Het_AB
Chr01	112623	Hom_B
Chr01	112626	Hom_B
Chr01	5278587	Hom_A
Chr01	12479439	Het_AB
Chr01	25279526	Hom_B
Chr01	35279704	Hom_A
Chr01	45280348	Het_AB
Chr02	1280730	Het_AB
Chr02	1680841	Hom_A
Chr02	2280919	Hom_B
Chr02	5281104	Hom_B
Chr02	15281539	Het_AB
Chr02	25281802	Het_AB
Chr02	35282172	Het_AB
Chr02	42282261	Hom_B`;

const defaultGeneData = 
`Satt300	Chr01	28572124	28572367	red
Satt429	Chr01	47217842	47218105	blue
Satt197	Chr02	8898878	8899050	green
Satt556	Chr02	38859467	38859623	yellow`;

// Data storage
let chrData = [];
let siteData = [];
let geneData = [];

// Function to download sample data
function downloadSampleData(dataType) {
    let content = '';
    let filename = '';
    
    switch(dataType) {
        case 'chromosome':
            content = defaultChrData;
            filename = 'sample_chromosome_data.txt';
            break;
        case 'site':
            content = defaultSiteData;
            filename = 'sample_site_data.txt';
            break;
        case 'gene':
            content = defaultGeneData;
            filename = 'sample_gene_data.txt';
            break;
        default:
            console.error('Unknown data type:', dataType);
            return;
    }
    
    // Create blob and download
    const blob = new Blob([content], { type: 'text/plain' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    window.URL.revokeObjectURL(url);
}

// Initialization
document.addEventListener('DOMContentLoaded', () => {
    // Initialize configuration parameters from DOM elements
    config.maxWidth = parseInt(document.getElementById('chrWidth').value);
    config.chromosomeHeight = parseInt(document.getElementById('chrHeight').value);
    config.chromosomeSpacing = parseInt(document.getElementById('chrSpacing').value);
    config.geneFontSize = parseInt(document.getElementById('geneFontSize').value);
    config.backgroundFontSize = parseInt(document.getElementById('backgroundFontSize').value);
    config.chrBorderWidth = parseInt(document.getElementById('chrBorderWidth').value);
    config.histogramBinSize = parseInt(document.getElementById('histogramBinSize').value);
    config.geneColorMode = document.getElementById('geneColorMode').value;
    config.geneAngle = parseInt(document.getElementById('geneAngle').value);
    config.alternateSides = document.getElementById('alternateSides').checked;
    
    // Initialize legend labels
    config.labels['Hom_A'] = document.getElementById('labelHomA').value;
    config.labels['Hom_B'] = document.getElementById('labelHomB').value;
    config.labels['Het_AB'] = document.getElementById('labelHetAB').value;
    config.labels['ParentA'] = document.getElementById('labelParentA').value;
    config.labels['ParentB'] = document.getElementById('labelParentB').value;
    
    // Bind file loading button
    document.getElementById('loadDataBtn').addEventListener('click', loadData);
    
    // Bind default data loading button
    document.getElementById('loadDefaultDataBtn').addEventListener('click', loadDefaultData);
    
    // Bind checkbox change events
    document.getElementById('showHomA').addEventListener('change', drawVisualization);
    document.getElementById('showHomB').addEventListener('change', drawVisualization);
    document.getElementById('showHetAB').addEventListener('change', drawVisualization);
    document.getElementById('showBackgroundAnalysis').addEventListener('change', drawVisualization);
    
    // Bind site display mode change event
    document.getElementById('siteDisplayMode').addEventListener('change', function() {
        config.displayMode = this.value;
        drawVisualization();
    });
    
    // Bind histogram window size selection event
    document.getElementById('histogramBinSize').addEventListener('change', function() {
        config.histogramBinSize = parseInt(this.value);
        drawVisualization();
    });
    
    // Bind chromosome size change events (changed to direct input)
    document.getElementById('chrWidth').addEventListener('change', function() {
        config.maxWidth = parseInt(this.value);
        drawVisualization();
    });
    // Add input event listener for real-time updates
    document.getElementById('chrWidth').addEventListener('input', function() {
        config.maxWidth = parseInt(this.value);
        drawVisualization();
    });
    
    document.getElementById('chrHeight').addEventListener('change', function() {
        config.chromosomeHeight = parseInt(this.value);
        drawVisualization();
    });
    // Add input event listener for real-time updates
    document.getElementById('chrHeight').addEventListener('input', function() {
        config.chromosomeHeight = parseInt(this.value);
        drawVisualization();
    });
    
    document.getElementById('chrSpacing').addEventListener('change', function() {
        config.chromosomeSpacing = parseInt(this.value);
        drawVisualization();
    });
    // Add input event listener for real-time updates
    document.getElementById('chrSpacing').addEventListener('input', function() {
        config.chromosomeSpacing = parseInt(this.value);
        drawVisualization();
    });
    
    // Bind font size change events
    document.getElementById('geneFontSize').addEventListener('change', function() {
        config.geneFontSize = parseInt(this.value);
        drawVisualization();
    });
    // Add input event listener for real-time updates
    document.getElementById('geneFontSize').addEventListener('input', function() {
        config.geneFontSize = parseInt(this.value);
        drawVisualization();
    });
    
    // Bind genetic background proportion font size change event
    document.getElementById('backgroundFontSize').addEventListener('change', function() {
        config.backgroundFontSize = parseInt(this.value);
        drawVisualization();
    });
    // Add input event listener for real-time updates
    document.getElementById('backgroundFontSize').addEventListener('input', function() {
        config.backgroundFontSize = parseInt(this.value);
        drawVisualization();
    });
    
    // Bind chromosome border width change event
    document.getElementById('chrBorderWidth').addEventListener('change', function() {
        config.chrBorderWidth = parseInt(this.value);
        drawVisualization();
    });
    
    // Bind chromosome color change event
    document.querySelectorAll('input[name="chrColor"]').forEach(radio => {
        radio.addEventListener('change', function() {
            // Immediately redraw visualization when color selection changes
            console.log('Chromosome color changed to:', this.value);
            drawVisualization();
        });
    });
    
    // Bind gene color mode selection event
    document.getElementById('geneColorMode').addEventListener('change', function() {
        config.geneColorMode = this.value;
        drawVisualization();
    });
    
    // Bind gene name display angle selection event
    document.getElementById('geneAngle').addEventListener('change', function() {
        config.geneAngle = parseInt(this.value);
        drawVisualization();
    });
    
    // Bind gene name alternating above/below display selection event
    document.getElementById('alternateSides').addEventListener('change', function() {
        config.alternateSides = this.checked;
        drawVisualization();
    });
    
    // Bind legend label change events
    document.getElementById('labelHomA').addEventListener('input', function() {
        config.labels['Hom_A'] = this.value;
        drawVisualization();
    });
    document.getElementById('labelHomB').addEventListener('input', function() {
        config.labels['Hom_B'] = this.value;
        drawVisualization();
    });
    document.getElementById('labelHetAB').addEventListener('input', function() {
        config.labels['Het_AB'] = this.value;
        drawVisualization();
    });
    
    // Bind parent label change events
    document.getElementById('labelParentA').addEventListener('input', function() {
        config.labels['ParentA'] = this.value;
        drawVisualization();
    });
    document.getElementById('labelParentB').addEventListener('input', function() {
        config.labels['ParentB'] = this.value;
        drawVisualization();
    });
    
    // Export-related settings have been updated to numeric input boxes, no longer need to update display values
    
    document.getElementById('exportBtn').addEventListener('click', exportVisualization);
});

// Load default data
async function loadDefaultData() {
    try {
        // Parse default data
        chrData = parseChrData(defaultChrData);
        siteData = parseSiteData(defaultSiteData);
        geneData = parseGeneData(defaultGeneData);
        
        // Simulate more data points for testing (generate 10,000 random sites per chromosome)
        const simulatedData = generateSimulatedData(chrData, 10000);
        siteData = [...siteData, ...simulatedData];
        
        // Draw visualization
        drawVisualization();
    } catch (error) {
        console.error('Error loading default data:', error);
        alert('Error loading default data: ' + error.message);
    }
}

// Generate simulated data for testing
function generateSimulatedData(chromosomes, pointsPerChromosome) {
    const simulatedData = [];
    const siteTypes = ['Hom_A', 'Hom_B', 'Het_AB'];
    
    chromosomes.forEach(chr => {
        for (let i = 0; i < pointsPerChromosome; i++) {
            // Randomly generate position and type
            const pos = Math.floor(Math.random() * chr.length);
            const typeIndex = Math.floor(Math.random() * 3);
            
            simulatedData.push({
                chrom: chr.chrom,
                pos: pos,
                siteType: siteTypes[typeIndex]
            });
        }
    });
    
    return simulatedData;
}

// Load user uploaded data
async function loadData() {
    try {
        // Get user uploaded files
        const chrFileInput = document.getElementById('chrFile');
        const siteFileInput = document.getElementById('siteFile');
        const geneFileInput = document.getElementById('geneFile');
        
        // Check if chromosome file is uploaded (chromosome file is required)
        if (!chrFileInput.files.length) {
            alert('Please upload at least the chromosome data file (chr.txt)');
            return;
        }
        
        // Reset data
        chrData = [];
        siteData = [];
        geneData = [];
        
        // Read chromosome data
        const chrText = await readFileContent(chrFileInput.files[0]);
        chrData = parseChrData(chrText);
        
        // Check if chromosome data is valid
        if (chrData.length === 0) {
            alert('Chromosome data is invalid or empty');
            return;
        }
        
        // Read site data (if available)
        if (siteFileInput.files.length) {
            const siteText = await readFileContent(siteFileInput.files[0]);
            siteData = parseSiteData(siteText);
        }
        
        // Read gene data (if available)
        if (geneFileInput.files.length) {
            const geneText = await readFileContent(geneFileInput.files[0]);
            geneData = parseGeneData(geneText);
        }
        
        // Draw visualization
        drawVisualization();
    } catch (error) {
        console.error('Error loading data:', error);
        alert('Error loading data: ' + error.message);
    }
}

// Read content from file
function readFileContent(file) {
    return new Promise((resolve, reject) => {
        const reader = new FileReader();
        reader.onload = (event) => resolve(event.target.result);
        reader.onerror = (error) => reject(error);
        reader.readAsText(file);
    });
}

// Parse chromosome data
function parseChrData(text) {
    const lines = text.trim().split('\n');
    return lines.map(line => {
        const [chrom, length] = line.split('\t');
        return { chrom, length: parseInt(length) };
    });
}

// Parse site data
function parseSiteData(text) {
    const lines = text.trim().split('\n');
    return lines.map(line => {
        const [chrom, pos, siteType] = line.split('\t');
        return { chrom, pos: parseInt(pos), siteType };
    });
}

// Parse gene data
function parseGeneData(text) {
    const lines = text.trim().split('\n');
    return lines.map(line => {
        const parts = line.split('\t');
        const gene = parts[0];
        const chr = parts[1];
        const startPos = parseInt(parts[2]);
        const endPos = parseInt(parts[3]);
        const color = parts.length >= 5 ? parts[4] : null; // Parse color if fifth column exists
        
        return { 
            gene, 
            chr, 
            startPos: startPos, 
            endPos: endPos,
            midPos: Math.floor((startPos + endPos) / 2),
            color: color // Store color information
        };
    });
}

// Calculate parental background percentage of chromosome
function calculateBackgroundPercentage(sites) {
    // Count different types of sites
    const counts = {
        'Hom_A': sites.filter(site => site.siteType === 'Hom_A').length,
        'Hom_B': sites.filter(site => site.siteType === 'Hom_B').length,
        'Het_AB': sites.filter(site => site.siteType === 'Het_AB').length
    };
    
    const totalSites = counts['Hom_A'] + counts['Hom_B'] + counts['Het_AB'];
    
    // Avoid division by zero
    if (totalSites === 0) {
        return { parentA: 0, parentB: 0 };
    }
    
    // Calculate parent A proportion: (2*n_A + n_AB)/(2*(n_A + n_B + n_AB))
    const parentA = (2 * counts['Hom_A'] + counts['Het_AB']) / (2 * totalSites);
    
    // Calculate parent B proportion: (2*n_B + n_AB)/(2*(n_A + n_B + n_AB))
    const parentB = (2 * counts['Hom_B'] + counts['Het_AB']) / (2 * totalSites);
    
    return { 
        parentA: parentA * 100, 
        parentB: parentB * 100,
        counts: counts
    };
}

// Determine gene color
function determineGeneColor(gene, sites) {
    // Determine returned color based on color mode
    
    // 1. Use custom color mode
    if (config.geneColorMode === 'custom' && gene.color) {
        return gene.color; // Directly use color specified in gene data
    }
    
    // 2. Use uniform black mode
    if (config.geneColorMode === 'black') {
        return '#333333'; // Uniformly use black
    }
    
    // 3. Use auto calculation mode (default)
    // Extract sites within gene range
    const geneSites = sites.filter(site => 
        site.pos >= gene.startPos && site.pos <= gene.endPos
    );
    
    // If no sites, return black
    if (geneSites.length === 0) {
        return '#333333';
    }
    
    // Check if contains Het_AB sites
    const hasHetAB = geneSites.some(site => site.siteType === 'Het_AB');
    if (hasHetAB) {
        return config.colors['Het_AB'];
    }
    
    // Check if only has Hom_A sites
    const hasHomA = geneSites.some(site => site.siteType === 'Hom_A');
    const hasHomB = geneSites.some(site => site.siteType === 'Hom_B');
    
    if (hasHomA && !hasHomB) {
        return config.colors['Hom_A'];
    }
    
    // Check if only has Hom_B sites
    if (hasHomB && !hasHomA) {
        return config.colors['Hom_B'];
    }
    
    // Default return black
    return '#333333';
}

// Draw visualization
function drawVisualization() {
    // Clear previous visualization
    d3.select('#visualization').html('');
    
    // If no chromosome data, do not draw
    if (chrData.length === 0) {
        return;
    }
    
    // Check if there is site data and gene data
    const hasSiteData = siteData.length > 0;
    const hasGeneData = geneData.length > 0;
    
    // Get site type filter status (only needed when there is site data)
    let filteredSiteData = [];
    if (hasSiteData) {
        const showHomA = document.getElementById('showHomA').checked;
        const showHomB = document.getElementById('showHomB').checked;
        const showHetAB = document.getElementById('showHetAB').checked;
        
        // Filter site data
        filteredSiteData = siteData.filter(site => {
            if (site.siteType === 'Hom_A' && !showHomA) return false;
            if (site.siteType === 'Hom_B' && !showHomB) return false;
            if (site.siteType === 'Het_AB' && !showHetAB) return false;
            return true;
        });
    }
    
    // Calculate maximum chromosome length
    const maxChrLength = Math.max(...chrData.map(d => d.length));
    
    // Dynamic calculation of left margin to ensure enough space for chromosome labels and gene names
    let dynamicLeftMargin = config.margin.left;
    
    // Consider chromosome label length
    const maxChrLabelLength = Math.max(...chrData.map(d => d.chrom.length));
    dynamicLeftMargin = Math.max(dynamicLeftMargin, maxChrLabelLength * 8 + 20);
    
    // If there is gene data, consider gene name length
    if (hasGeneData) {
        const maxGeneLabelLength = Math.max(...geneData.map(d => d.gene.length));
        dynamicLeftMargin = Math.max(dynamicLeftMargin, maxGeneLabelLength * 6 + 40);
    }
    
    // Ensure left margin is not less than minimum value
    dynamicLeftMargin = Math.max(dynamicLeftMargin, 120);
    
    // Calculate visualization dimensions
    const width = Math.min(config.maxWidth, window.innerWidth - 40);
    const height = chrData.length * config.chromosomeSpacing + config.margin.top + config.margin.bottom;
    
    // Create SVG element, add viewBox to ensure full display
    const svg = d3.select('#visualization')
        .append('svg')
        .attr('width', width)
        .attr('height', height)
        .attr('viewBox', `0 0 ${width} ${height}`)
        .attr('preserveAspectRatio', 'xMinYMin meet')
        .attr('id', 'chromosomeVisualization');
    
    // Create main drawing area, use dynamic left margin
    const g = svg.append('g')
        .attr('transform', `translate(${dynamicLeftMargin},${config.margin.top})`);
    
    // Create x scale, use dynamic left margin
    const xScale = d3.scaleLinear()
        .domain([0, maxChrLength])
        .range([0, width - dynamicLeftMargin - config.margin.right]);
    
    // For each chromosome, create a group
    const chromosomes = g.selectAll('.chromosome-group')
        .data(chrData)
        .enter()
        .append('g')
        .attr('class', 'chromosome-group')
        .attr('transform', (d, i) => `translate(0,${i * config.chromosomeSpacing})`);
    
    // Get selected chromosome color
    const chrColor = document.querySelector('input[name="chrColor"]:checked').value;
    
    // Define clip path to ensure histogram only displays within chromosome
    chromosomes.each(function(d, i) {
        const clipPathId = `clip-chromosome-${i}`;
        
        d3.select(this).append('clipPath')
            .attr('id', clipPathId)
            .append('rect')
            .attr('x', 0)
            .attr('y', 0)
            .attr('width', xScale(d.length))
            .attr('height', config.chromosomeHeight)
            .attr('rx', config.chromosomeHeight / 2)
            .attr('ry', config.chromosomeHeight / 2);
    });
    
    // Draw chromosome fill part
    chromosomes.append('rect')
        .attr('class', 'chromosome-fill')
        .attr('x', 0)
        .attr('y', 0)
        .attr('width', d => xScale(d.length))
        .attr('height', config.chromosomeHeight)
        .attr('rx', config.chromosomeHeight / 2)
        .attr('ry', config.chromosomeHeight / 2)
        .attr('fill', chrColor); // Use selected chromosome color
    
    // If there is site data, draw histogram after chromosome fill
    if (hasSiteData) {
        chrData.forEach((chr, chrIndex) => {
            const chrSites = filteredSiteData.filter(site => site.chrom === chr.chrom);
            const chromosomeGroup = chromosomes.filter((d, i) => i === chrIndex);
            
            // Create histogram group, apply clip path
            const histogramGroup = chromosomeGroup.append('g')
                .attr('class', 'histogram-group')
                .attr('clip-path', `url(#clip-chromosome-${chrIndex})`);
            
            // Draw histogram
            drawHistogram(histogramGroup, chrSites, chr.length, xScale);
        });
    }
    
    // Draw chromosome border (draw after histogram to ensure border is above)
    chromosomes.append('rect')
        .attr('class', 'chromosome')
        .attr('x', 0)
        .attr('y', 0)
        .attr('width', d => xScale(d.length))
        .attr('height', config.chromosomeHeight)
        .attr('rx', config.chromosomeHeight / 2)
        .attr('ry', config.chromosomeHeight / 2)
        .attr('fill', 'none') // Border does not fill color
        .attr('stroke', '#000000') // Add border color
        .attr('stroke-width', config.chrBorderWidth) // Use custom border width
    
    // Add chromosome label to ensure enough left margin
    chromosomes.append('text')
        .attr('class', 'chromosome-label')
        .attr('x', -15) // Add offset to ensure not truncated
        .attr('y', config.chromosomeHeight / 2)
        .attr('dy', '0.35em')
        .attr('text-anchor', 'end')
        .text(d => d.chrom);
    
    // If there is site data and show parental background proportion analysis
    if (hasSiteData) {
        // Get checkbox status of whether to show background proportion analysis
        const showBackgroundAnalysis = document.getElementById('showBackgroundAnalysis').checked;
        
        chrData.forEach((chr, chrIndex) => {
            const chrSites = filteredSiteData.filter(site => site.chrom === chr.chrom);
            const percentages = calculateBackgroundPercentage(chrSites);
            
            const chromosomeGroup = chromosomes.filter((d, i) => i === chrIndex);
            
            // If choose to show background proportion analysis, add parental A and parental B proportions
            if (showBackgroundAnalysis) {
                // Add parental A proportion
                chromosomeGroup.append('text')
                    .attr('class', 'background-percentage')
                    .attr('x', xScale(chr.length) + 10) // Reduce offset to let text closer to chromosome
                    .attr('y', config.chromosomeHeight / 2 - 3) // Move down a little
                    .attr('text-anchor', 'start')
                    .attr('fill', config.colors['Hom_A'])
                    .attr('font-size', `${config.backgroundFontSize}px`) // Use configured font size
                    .text(`${config.labels['ParentA']} ${percentages.parentA.toFixed(2)}%`);
                
                // Add parental B proportion
                chromosomeGroup.append('text')
                    .attr('class', 'background-percentage')
                    .attr('x', xScale(chr.length) + 10) // Reduce offset to let text closer to chromosome
                    .attr('y', config.chromosomeHeight / 2 + 7) // Move down a little
                    .attr('text-anchor', 'start')
                    .attr('fill', config.colors['Hom_B'])
                    .attr('font-size', `${config.backgroundFontSize}px`) // Use configured font size
                    .text(`${config.labels['ParentB']} ${percentages.parentB.toFixed(2)}%`);
            }
        });
    }
    
    // If there is gene data, draw gene markers
    if (hasGeneData) {
        chrData.forEach((chr, chrIndex) => {
            const chrSites = hasSiteData ? filteredSiteData.filter(site => site.chrom === chr.chrom) : [];
            const chrGenes = geneData.filter(gene => gene.chr === chr.chrom);
            
            // Calculate each gene's position above or below for alternating display mode
            if (config.alternateSides) {
                chrGenes.forEach((gene, idx) => {
                    gene.position = idx % 2 === 0 ? 'top' : 'bottom';
                });
            }
            
            const geneGroup = chromosomes.filter((d, i) => i === chrIndex)
                .selectAll('.gene-group')
                .data(chrGenes)
                .enter()
                .append('g')
                .attr('class', 'gene-group');
            
            // Draw connection line
            geneGroup.append('line')
                .attr('class', 'gene-line')
                .attr('x1', d => xScale(d.midPos))
                .attr('y1', d => config.alternateSides && d.position === 'bottom' ? 
                      config.chromosomeHeight + 5 : -5) // Determine start point based on display mode
                .attr('x2', d => xScale(d.midPos))
                .attr('y2', d => config.alternateSides && d.position === 'bottom' ? 
                      config.chromosomeHeight + config.geneConnectorOffset : 
                      -config.geneConnectorOffset) // Determine end point based on display mode
                .attr('stroke', d => determineGeneColor(d, chrSites)); // Use color determination function
            
            // Add gene label to ensure not exceed left boundary
            geneGroup.append('text')
                .attr('class', 'gene-label')
                .attr('x', d => {
                    const baseX = xScale(d.midPos) + 1;
                    // Ensure gene label does not exceed left boundary
                    return Math.max(baseX, 5);
                })
                .attr('y', d => config.alternateSides && d.position === 'bottom' ? 
                      config.chromosomeHeight + 12 : -6) // Determine y position based on display mode
                .attr('text-anchor', d => {
                    const baseX = xScale(d.midPos) + 1;
                    // If position adjusted to left boundary, change alignment
                    return baseX < 5 ? 'start' : 'start';
                })
                .attr('fill', d => determineGeneColor(d, chrSites))
                .attr('font-size', `${config.geneFontSize}px`)
                .attr('transform', d => {
                    // Calculate rotation angle
                    const baseX = xScale(d.midPos) + 1;
                    const x = Math.max(baseX, 5);
                    const y = config.alternateSides && d.position === 'bottom' ? 
                              config.chromosomeHeight + 12 : -6;
                    
                    // Apply user-specified angle, use opposite angle for bottom part
                    let angle = config.geneAngle;
                    if (config.alternateSides && d.position === 'bottom') {
                        angle = -angle; // Bottom label uses opposite angle
                    }
                    
                    return `rotate(${angle}, ${x}, ${y})`;
                })
                .text(d => d.gene);
                
            // Add mouseover event for gene group
            geneGroup.on('mouseover', function() {
                d3.select(this).select('.gene-label').style('font-weight', '600');
            }).on('mouseout', function() {
                d3.select(this).select('.gene-label').style('font-weight', '500');
            });
        });
    }
    
    // Draw beautified horizontal axis, use dynamic left margin
    drawChromosomeAxis(svg, xScale, width, height, maxChrLength, dynamicLeftMargin);
    
    // If there is site data, add legend
    if (hasSiteData) {
        const legend = svg.append('g')
            .attr('class', 'legend')
            .attr('transform', `translate(${width / 2}, ${height - 15})`);
        
        const legendItems = [
            { type: 'Hom_A', label: config.labels['Hom_A'] },
            { type: 'Hom_B', label: config.labels['Hom_B'] },
            { type: 'Het_AB', label: config.labels['Het_AB'] }
        ];
        
        const legendG = legend.selectAll('.legend-item')
            .data(legendItems)
            .enter()
            .append('g')
            .attr('class', 'legend-item')
            .attr('transform', (d, i) => `translate(${i * 90 - 135}, 0)`);
        
        legendG.append('circle')
            .attr('r', 5)
            .attr('fill', d => config.colors[d.type]);
        
        legendG.append('text')
            .attr('x', 10)
            .attr('y', 4)
            .text(d => d.label);
    }
}

// Draw chromosome-specific horizontal axis
function drawChromosomeAxis(svg, xScale, width, height, maxChrLength, dynamicLeftMargin = config.margin.left) {
    const axisGroup = svg.append('g')
        .attr('transform', `translate(${dynamicLeftMargin},${height - config.margin.bottom})`);
    
    // Main axis line
    axisGroup.append('line')
        .attr('class', 'axis-main-line')
        .attr('x1', 0)
        .attr('y1', 0)
        .attr('x2', xScale(maxChrLength) + 20) // Extend axis line to display complete arrow
        .attr('y2', 0)
        .attr('stroke', '#333')
        .attr('stroke-width', 2); // Thick main axis line
    
    // Automatically calculate suitable tick interval
    function calculateTickInterval(maxLength) {
        const maxTicks = 10; // Maximum display 10 ticks
        const lengthInMb = maxLength / 1000000;
        
        // Optional interval values (Mb)
        const intervals = [1, 2, 5, 10, 20, 25, 50, 100, 200, 500, 1000];
        
        // Find suitable interval to ensure tick count does not exceed maxTicks
        for (let interval of intervals) {
            if (lengthInMb / interval <= maxTicks) {
                return interval * 1000000; // Convert to bp
            }
        }
        
        // If none suitable, use maximum interval
        return intervals[intervals.length - 1] * 1000000;
    }
    
    const majorTickStep = calculateTickInterval(maxChrLength);
    const majorTickCount = Math.ceil(maxChrLength / majorTickStep);
    
    // Generate main ticks and labels
    for (let i = 0; i <= majorTickCount; i++) {
        const pos = i * majorTickStep;
        if (pos > maxChrLength) break; // Ensure not exceed longest chromosome length
        
        const xPos = xScale(pos);
        
        // Draw main tick line
        axisGroup.append('line')
            .attr('class', 'axis-major-tick')
            .attr('x1', xPos)
            .attr('y1', 0)
            .attr('x2', xPos)
            .attr('y2', 8)
            .attr('stroke', '#333')
            .attr('stroke-width', 1.5); // Thick tick line
        
        // Add text label, select appropriate unit based on size
        const labelText = pos >= 1000000 ? `${pos / 1000000}Mb` : `${pos / 1000}kb`;
        axisGroup.append('text')
            .attr('class', 'axis-major-label')
            .attr('x', xPos)
            .attr('y', 25)
            .attr('text-anchor', 'middle')
            .attr('font-size', '12px')
            .attr('fill', '#333')
            .text(labelText);
    }
    
    // Add complete arrow
    axisGroup.append('path')
        .attr('class', 'axis-arrow')
        .attr('d', `M${xScale(maxChrLength)},0 L${xScale(maxChrLength) + 15},0 L${xScale(maxChrLength) + 15},-5 L${xScale(maxChrLength) + 25},0 L${xScale(maxChrLength) + 15},5 L${xScale(maxChrLength) + 15},0 Z`)
        .attr('fill', '#333');
}

// Draw histogram
function drawHistogram(histogramGroup, sites, chrLength, xScale) {
    // Calculate histogram bar
    const binWidth = config.histogramBinSize;
    const binCount = Math.ceil(chrLength / binWidth);
    
    // Group sites by site type
    const sitesByType = {
        'Hom_A': sites.filter(s => s.siteType === 'Hom_A'),
        'Hom_B': sites.filter(s => s.siteType === 'Hom_B'),
        'Het_AB': sites.filter(s => s.siteType === 'Het_AB')
    };
    
    // Calculate site count in each bar
    const histogramData = [];
    
    for (let i = 0; i < binCount; i++) {
        const startPos = i * binWidth;
        const endPos = Math.min((i + 1) * binWidth, chrLength);
        
        // Calculate site count in each type in this interval
        const countHomA = sitesByType['Hom_A'].filter(s => s.pos >= startPos && s.pos < endPos).length;
        const countHomB = sitesByType['Hom_B'].filter(s => s.pos >= startPos && s.pos < endPos).length;
        const countHetAB = sitesByType['Het_AB'].filter(s => s.pos >= startPos && s.pos < endPos).length;
        
        if (countHomA > 0) {
            histogramData.push({
                startPos,
                endPos,
                count: countHomA,
                type: 'Hom_A'
            });
        }
        
        if (countHomB > 0) {
            histogramData.push({
                startPos,
                endPos,
                count: countHomB,
                type: 'Hom_B'
            });
        }
        
        if (countHetAB > 0) {
            histogramData.push({
                startPos,
                endPos,
                count: countHetAB,
                type: 'Het_AB'
            });
        }
    }
    
    // Calculate maximum bar height
    const maxCount = Math.max(...histogramData.map(d => d.count), 1); // Avoid division by 0
    const barMaxHeight = config.chromosomeHeight * 0.8; // Slightly reduce height to ensure not exceed chromosome edge
    
    // Draw histogram bar
    histogramGroup.selectAll('.histogram-bar')
        .data(histogramData)
        .enter()
        .append('rect')
        .attr('class', d => `histogram-bar histogram-bar-${d.type}`)
        .attr('x', d => xScale(d.startPos))
        .attr('y', d => {
            // Calculate bar height
            const barHeight = d.count / maxCount * barMaxHeight;
            return config.chromosomeHeight / 2 - barHeight / 2;
        })
        .attr('width', d => Math.max(1, xScale(d.endPos) - xScale(d.startPos) - 0.5)) // Leave a little space
        .attr('height', d => d.count / maxCount * barMaxHeight)
        .attr('fill', d => config.colors[d.type])
        .attr('opacity', 0.8)
        .attr('rx', 1) // Give bar rounded corners
        .attr('ry', 1);
}

// Export visualization image
function exportVisualization() {
    try {
        const format = document.getElementById('exportFormat').value;
        const quality = parseFloat(document.getElementById('exportQuality').value);
        const scale = parseFloat(document.getElementById('exportScale').value);
        
        // Get selected chromosome color
        const chrColor = document.querySelector('input[name="chrColor"]:checked').value;
        
        const visualizationElement = document.getElementById('chromosomeVisualization');
        
        if (!visualizationElement) {
            alert('Please load data first to generate visualization image');
            return;
        }
        
        // Create cloned version of SVG for export, avoid modifying original SVG
        const clonedSvg = visualizationElement.cloneNode(true);
        
        // Ensure SVG has correct viewBox to prevent content truncation
        const svgWidth = visualizationElement.getAttribute('width');
        const svgHeight = visualizationElement.getAttribute('height');
        
        // If original SVG does not have viewBox, add one
        if (!clonedSvg.getAttribute('viewBox')) {
            clonedSvg.setAttribute('viewBox', `0 0 ${svgWidth} ${svgHeight}`);
        }
        
        // Ensure preserveAspectRatio is correctly set
        clonedSvg.setAttribute('preserveAspectRatio', 'xMinYMin meet');
        
        // Chromosome fill color has been set during drawing, no need to set again
        // Border style has been set in drawVisualization, no need to repeat
        
        // Export method
        if (format === 'svg') {
            // SVG export
            const serializer = new XMLSerializer();
            const svgString = serializer.serializeToString(clonedSvg);
            const svgBlob = new Blob([svgString], {type: 'image/svg+xml'});
            downloadBlob(svgBlob, 'Chromosome Visualization.svg');
        } else {
            const canvas = document.createElement('canvas');
            const serializer = new XMLSerializer();
            const svgString = serializer.serializeToString(clonedSvg);
            const svgBlob = new Blob([svgString], {type: 'image/svg+xml;charset=utf-8'});
            const DOMURL = window.URL || window.webkitURL || window;
            const url = DOMURL.createObjectURL(svgBlob);
            
            const img = new Image();
            img.onload = function() {
                canvas.width = visualizationElement.width.baseVal.value * scale;
                canvas.height = visualizationElement.height.baseVal.value * scale;
                const ctx = canvas.getContext('2d');
                
                // If JPEG format, need to fill white background first
                if (format === 'jpeg') {
                    ctx.fillStyle = '#ffffff';
                    ctx.fillRect(0, 0, canvas.width, canvas.height);
                }
                
                ctx.scale(scale, scale);
                ctx.drawImage(img, 0, 0);
                DOMURL.revokeObjectURL(url);
                
                let imgType = 'image/png';
                let fileName = 'Chromosome Visualization.png';
                
                if (format === 'jpeg') {
                    imgType = 'image/jpeg';
                    fileName = 'Chromosome Visualization.jpeg';
                }
                
                canvas.toBlob(function(blob) {
                    downloadBlob(blob, fileName);
                }, imgType, quality);
            };
            
            img.onerror = function(error) {
                console.error('Image loading error:', error);
                alert('Error exporting image: Unable to load SVG as image');
            };
            
            img.src = url;
        }
    } catch (error) {
        console.error('Error during export:', error);
        alert('Error during export: ' + (error.message || 'Unknown error'));
    }
}

// Download Blob
function downloadBlob(blob, fileName) {
    const link = document.createElement('a');
    link.href = URL.createObjectURL(blob);
    link.download = fileName;
    document.body.appendChild(link);
    link.click();
    setTimeout(function() {
        document.body.removeChild(link);
        URL.revokeObjectURL(link.href);
    }, 100);
}