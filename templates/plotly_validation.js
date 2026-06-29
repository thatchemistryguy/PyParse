function refreshValidation() {
    let min_time = 100
    let max_time = 0
    
    validation_plot = []
    let max_area_green = 0
    
    if(cpData[this_cmp_no]["final_result"]["green"].length > 0) {
        max_area_green = peaks[cpData[this_cmp_no]["final_result"]["green"].reduce((prev, current) => (prev && peaks[prev]["area"] > peaks[current]["area"]) ? prev : current)]["area"]
    }

    let max_area_discarded = 0
    if(cpData[this_cmp_no]["final_result"]["discarded"].length > 0) {
        max_area_discarded = peaks[cpData[this_cmp_no]["final_result"]["discarded"].reduce((prev, current) => (prev && peaks[prev]["area"] > peaks[current]["area"]) ? prev : current)]["area"]
    }

    let max_area_discarded_cluster = 0
    cpData[this_cmp_no]["discarded_clusters"].forEach(function(j, index) {
        for (const[key, value] of Object.entries(j)) {
            if(value.length > 0) {
                let max_in_this_round = value.reduce((prev, current) => (prev && peaks[prev]["area"] > peaks[current]["area"]) ? prev : current)
                if(peaks[max_in_this_round]["area"] > max_area_discarded_cluster) {
                    max_area_discarded_cluster = peaks[max_in_this_round]["area"]
                }
            }
        }
    })

    

    let max_area_global = Math.max(max_area_green, max_area_discarded, max_area_discarded_cluster)


    cpData[this_cmp_no]["final_result"]["green"].forEach(function(i, index) {
        let x = peaks[i]["well"]
        let y = peaks[i]["time"]
        let marker_size = 15 * (peaks[i]["area"] / max_area_global)
        let marker_colour = "blue"
        let marker_shape = "circle"
        
        let new_point = new new_validation_point(x, y, marker_size, marker_colour, marker_shape)
        validation_plot.push({...new_point})

        if(y < min_time){
            min_time = y
        }
        if(y > max_time){
            max_time = y
        }
    })

    cpData[this_cmp_no]["final_result"]["discarded"].forEach(function(j, index) {
        let x = peaks[j]["well"]
        let y = peaks[j]["time"]
        let marker_size = 15 * (peaks[j]["area"] / max_area_global)
        let marker_colour = "yellow"
        let marker_shape = "circle"

        let new_point = new new_validation_point(x, y, marker_size, marker_colour, marker_shape)
        validation_plot.push({...new_point})

        if(y < min_time){
            min_time = y
        }
        if(y > max_time){
            max_time = y
        }
    })

    let marker_symbols = ["triangle-up", "triangle-down", "triangle-ne", "triangle-se", "triangle-sw", "triangle-nw", "pentagon", "hexagon", "hexagon2", "octagon"]
    let counter = -1
    cpData[this_cmp_no]["discarded_clusters"].forEach(function(j, index) {
        counter++
        if(counter == marker_symbols.length) {
            counter = 0
        }
        for (const[key, k] of Object.entries(j)) {
            k.forEach(function(peakID, bindex) {
                let x = peaks[peakID]["well"]
                let y = peaks[peakID]["time"]
                let marker_size = 15 * (peaks[peakID]["area"] / max_area_global)
                let marker_colour = "red"
                let marker_shape = marker_symbols[counter]

                let new_point = new new_validation_point(x, y, marker_size, marker_colour, marker_shape)
                validation_plot.push({...new_point})

                if(y < min_time){
                    min_time = y
                }
                if(y > max_time){
                    max_time = y
                }
            })
        }
        
    })


    min_time = min_time - 0.05
    max_time = max_time + 0.05
    
    validation_layout["yaxis"]["range"] = [min_time, max_time]
    Plotly.newPlot("validation-plot-space", validation_plot, validation_layout)
    setValidationClickFns()

}

class new_validation_point {
    constructor(x, y, marker_size, marker_colour, marker_shape) {
        this.type = "scatter"
        this.mode = "markers"
        this.marker = {
            color: marker_colour,
            size: marker_size,
            symbol: marker_shape
        }
        this.x = [x]
        this.y = [y]
        this.showlegend = false
    }
}



let validation_layout = {
    xaxis: {
        zeroline: false,
        range: [0, 1],
        rangemode: "tozero",
        title: {
            text: "Well Number",
            font: {
                size: 14,
            }
        }
    },
    yaxis: {
        zeroline: false,
        title: {
            text: 'Time /min',
            font: {
                size: 14,
            }
        }
    },
    height:800,
    autosize: true,
}
let validation_plot = []

//On clicking a point on the validation plot, change the chroma_plot to show that particular well. 
function setValidationClickFns() {
    let myPlot = document.getElementById('validation-plot-space')
    myPlot.on('plotly_click', function(clickdata){
        if(clickdata.points.length > 0 && clickdata.points[0].y !== undefined) {
            let filtered_samples = samples.filter((x) => x["well"] == clickdata.points[0].x)
            this_sample_name = filtered_samples[0]["filename"]
            this_sample_no = samples.findIndex(x => x["filename"] == this_sample_name)

            //Find the nearest peakID to this retention time
            let filtered_lc = lcData.filter((x) => x["filename"] == this_sample_name)
            let closest_peak = filtered_lc.reduce((prev, current) => (prev && Math.abs(prev["time"]-clickdata.points[0].y) < Math.abs(current["time"]-clickdata.points[0].y)) ? prev : current, {"time": 1000})
            refreshChroma(closest_peak["peakID"])


            //Set the name of the sample currently shown. 
            $(' #curr_sample_name ').html(cpData[this_cmp_no]["name"] + '; well: ' + samples[this_sample_no]["well_readable"] + '.')
        }
    });
}

