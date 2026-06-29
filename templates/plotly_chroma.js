
let plotting = {
    name: "",
    x: [],
    y: [],
    type: "scatter",
    mode: 'lines+text',
    line: {
        color: '#0000FF',
        width: 1,
        shape: 'spline', 
    },
    xaxis: "x1",
    yaxis: "y1",
    customdata: [],
    hoveron: "points+fills",
    hovertemplate: '%{customdata}',
    showlegend: false,
}

let overlay = {
    name: "",
    x: [],
    y: [],
    type: "scatter",
    fill: "tonexty",
    fillcolor:"#0000FF",
    mode: 'lines+text',
    line: {
        color: '#0000FF',
        width: 1,
        shape: 'spline', 
    },
    xaxis: "x1",
    yaxis: "y1",
    customdata: [],
    hoveron: "points+fills",
    hovertemplate: '%{customdata}',
    showlegend: false

}

let uv_plotting = {
    name: "",
    x: [],
    y: [],
    type: "scatter", 
    mode: 'lines+text',
    line: {
        color: '#0000FF',
        width: 1,
        shape: 'spline', 
    },
    xaxis: "x2",
    yaxis: "y2",
    showlegend: false
}

let ms_plus = {
    name: "",
    x: [],
    y: [],
    type: "bar", 
    marker: {
        color: 'rgba(50,171, 96, 1)',
    },
    width: 1,
    xaxis: "x3",
    yaxis: "y3",
    showlegend: false
}

let ms_minus = {
    name: "",
    x: [],
    y: [],
    type: "bar", 
    marker: {
        color: 'rgba(50,171, 96, 1)',
    },
    width: 1,
    xaxis: "x4",
    yaxis: "y4",
    showlegend: false
}

let buttons_to_remove = ["pan2d", "select2d", "lasso2d", "zoomIn2d", "zoomOut2d", "autoScale2d", "hoverClosestGl2d", "hoverClosestPie", "toggleHover", "toImage", "sendDataToCloud", "toggleSpikelines", "resetViewMapbox"]
let layout = {
    xaxis: {
        zeroline: false,
        domain: [0, 1],
        anchor: "y1", 
        rangemode: "tozero",
        title: {
            text: 'Time /min',
            font: {
                size: 14,
            }
        }
    },
    yaxis: {
        zeroline: false,
        domain: [0, 0.55],
        anchor: "x1",
        title: {
            text: 'DAD Chromatogram',
            font: {
                size: 14,
            }
        }
    },
    xaxis2: {
        zeroline: false,
        domain: [0, 1],
        anchor: "y2",
    },
    yaxis2: {
        domain: [0.6, 0.7],
        anchor: "x2",
        title: {
            text: 'UV',
            font: {
                size: 14,
            }
        }
    },
    xaxis3: {
        zeroline: false,
        domain: [0, 1],
        anchor: "y3",
        range: [0,1000],
    },
    yaxis3: {
        domain: [0.75, 0.85],
        anchor: "x3",
        title: {
            text: 'MS+',
            font: {
                size: 14,
            }
        }
    },
    xaxis4: {
        zeroline: false,
        domain: [0, 1],
        anchor: "y4",
        range: [0, 1000],
    },
    yaxis4: {
        domain: [0.9, 1],
        anchor: "x4",
        title: {
            text: 'MS-',
            font: {
                size: 14,
            }
        }
    },
    grid: {rows: 3, columns: 1, pattern:"independent"},
    height:800,
    autosize: true,
}

function setClickFns() {
    let myPlot = document.getElementById('chroma-plot-space')
    myPlot.on('plotly_click', function(clickdata){
        if(clickdata.points.length > 0 && clickdata.points[0].pointNumber !== undefined) {
            let peakID = data_text[clickdata.points[0].pointNumber]
            if(peakID != -1) {
                refreshChroma(peakID)
            }
        }
        else if(clickdata.points[0].pointNumber === undefined) {
            let cur_x = clickdata.event.x
            let min_x = clickdata.points[0].xaxis.range[0]
            let max_x = clickdata.points[0].xaxis.range[1]
            
            //Upon click, the event target is the rect element which 
            //contains just the plotting area (doesn't include the axis)
            let svg_x = clickdata.event.target.x.animVal.value
            let svg_width = clickdata.event.target.width.animVal.value

            //Get the position of the parent svg element that contains all of the plotly drawing
            let parent_x = $(clickdata.event.srcElement.viewportElement).position().left
            
            //Find how far along the x-axis the click occurred. 
            let distance_from_left = (cur_x - svg_x - parent_x) / svg_width

            //graph_x is the x-value on the graphed at the point the click occurred. 
            let graph_x = min_x + distance_from_left * (max_x - min_x)
            
            //Find the nearest peakID to this position
            let filtered_lc = lcData.filter((x) => x["filename"] == this_sample_name)
            let closest_peak = filtered_lc.reduce((prev, current) => (prev && Math.abs(prev["time"]-graph_x) < Math.abs(current["time"]-graph_x)) ? prev : current, {"time": 1000})
            refreshChroma(closest_peak["peakID"])
        }
        
    });
}

function refreshChroma(peakID) {
    let times = []
    let heights = []
    let annotations = []
    let heights_amended = []
    let ms_minusx = []
    let ms_minusy = []
    
    let ms_plusx = []
    let ms_plusy = []

    let uv_x = []
    let uv_y = []
    data_text = []
    let custom_hover_text = []



    let filtered_lc = lcData.filter((value) => value["filename"] == this_sample_name)
    let filtered_ms = msData.filter((value) => value["filename"] == this_sample_name && value["peakID"] == peakID)
    let peak_data = peaks.filter((value) => value["filename"] == this_sample_name && value["peakID"] == peakID)
    let filtered_uv = uvData.filter((value) => value["filename"] == this_sample_name && value["peakID"] == peakID)
    let filtered_sample = samples.filter((value) => value["filename"] == this_sample_name)

    if(peak_data.length > 0) {
        peak_data = peak_data[0]
        
        let start_x = 0
        let start_y = 0
        $.each(filtered_lc, function(index, value) {
            times.push(value["time"])
            heights.push(value["height"])
            data_text.push(value["peakID"])
            if(value["peakID"] != -1) {
                custom_hover_text.push("Peak " + value["peakID"] + " <br>" + value["time"] + " min. ")
                
                
            }
            else {
                custom_hover_text.push("No peak here. <br>" + value["time"] + " min. ")
                //if this point isn't currently labelled a peak, set this as the new start_x for
                //when there is a peak
                start_x = value["time"]
                start_y = value["height"]
            }
            
            if(value["time"] <= peak_data["pEnd"] && value["time"] >= peak_data["pStart"]){
                //Get height at start of peak

                //Get height of the trace the next time there is a point
                //that isn't part of a peak. 
                let filter2 = filtered_lc.filter(function(bvalue, bindex) {
                    if(bindex > index && bvalue["peakID"] == -1) {
                        return true
                    }
                })
                let end_x = 0
                let end_y = 0
                if(filter2.length > 0) {
                    end_x = filter2[0]["time"]
                    end_y = filter2[0]["height"]
                }
                else {
                    end_x = filtered_lc[filtered_lc.length-1]["time"]
                    end_y = filtered_lc[filtered_lc.length-1]["height"]
                }
                let new_height = start_y + ((value["time"] - start_x)/(end_x - start_x) * (end_y - start_y))
                heights_amended.push(new_height)
            }
            else {
                heights_amended.push(value["height"])
            }

        })

        $.each(peaks.filter((x) => x["filename"] == this_sample_name), function(index, value){
            
            //Get height of the peak at its maxima
            let diffs = filtered_lc.map((x) => Math.abs(x["time"] - value["time"]))

            let x_index = diffs.indexOf(Math.min(...diffs))

            annotations.push({
                x: value["time"],
                y: heights[x_index] + 20,
                text: '(' + value["peakID"] + ')<br>' + value["time"] + '<br>' + value["area"] + '%',
                xref: "x",
                yref: "y",
                showarrow:false,
                font: {size: 10}
            })
        })

        //Set the x-axis to match the length of the run
        layout["xaxis"]["range"] = [0, Math.max(...filtered_lc.map((x) => x["time"]))]

        //Find max value of ms data
        let max_ms = filtered_ms.map((x) => x["MSvalue"])
        max_ms.push(650)

        layout["xaxis3"]["range"] = [0, Math.max(...max_ms)]
        layout["xaxis4"]["range"] = [0, Math.max(...max_ms)]

        layout["annotations"] = annotations

        plotting["x"] = times
        plotting["y"] = heights
        plotting["customdata"] = custom_hover_text

        overlay["x"] = times
        overlay["y"] = heights_amended
        overlay["customdata"] = custom_hover_text


        
        $.each(filtered_ms, function(index, value) {
            if(value["MStype"] == "-") {
                ms_minusx.push(value["MSvalue"])
                ms_minusy.push(value["MSintensity"])
                
                if(parseFloat(value["MSintensity"]) > 50) {
                    annotations.push({
                        x: parseFloat(value["MSvalue"]),
                        y: 100,
                        text: value["MSvalue"],
                        xref: "x4",
                        yref: "y4",
                        showarrow:true,
                        font: {size: 10},
                        ax: 0,
                        ay: -5
                    })
                }
            }
            else {
                ms_plusx.push(value["MSvalue"])
                ms_plusy.push(value["MSintensity"])
                
                if(parseFloat(value["MSintensity"]) > 50) {
                    annotations.push({
                        x: parseFloat(value["MSvalue"]),
                        y: 100,
                        text: value["MSvalue"],
                        xref: "x3",
                        yref: "y3",
                        showarrow:true,
                        font: {size: 10},
                        ax: 0,
                        ay: -5
                    })
                }
            }
        })

        ms_minus["x"] = ms_minusx
        ms_minus["y"] = ms_minusy

        ms_plus["x"] = ms_plusx
        ms_plus["y"] = ms_plusy
        
        $.each(filtered_uv, function(index, value) {
            uv_x.push(value["UVvalue"])
            uv_y.push(value["UVintensity"])

        })

        uv_plotting["x"] = uv_x
        uv_plotting["y"] = uv_y

        Plotly.redraw("chroma-plot-space")
        $(' #titleh1 ').html(this_sample_name + "; well: " + filtered_sample[0]["well"])
        $(' #methodh2 ').html(filtered_sample[0]["methodID"])
    }
}
