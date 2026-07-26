function changeHeatmap() {
    heatmap_plot["x"] = heatmapData[this_heatmap]["x_values"]
    heatmap_plot["y"] = heatmapData[this_heatmap]["y_values"]
    heatmap_plot["z"] = heatmapData[this_heatmap]["z_values"]
    let max_value = 0
    heatmapData[this_heatmap]["z_values"].forEach(function(row, row_index) {
        row.forEach(function(i, cell_index) {
            if(i>max_value) {
                max_value = i
            }
        })
    })

    heatmap_layout["annotations"] = []
    heatmapData[this_heatmap]["z_values"].forEach(function(row, row_index) {
        row.forEach(function(i, cell_index) {
            let textColour = ""

            if(i/max_value > 0.4) {
                textColour = "black"
            }
            else {
                textColour = "white"
            }
            let result = {
                xref: 'x1',
                yref: 'y1',
                x: heatmapData[this_heatmap]["x_values"][cell_index],
                y: heatmapData[this_heatmap]["y_values"][row_index],
                text: i,
                font: {
                    family: 'Arial',
                    size: 12,
                    color: textColour
                },
                showarrow: false,
            }
            heatmap_layout["annotations"].push(result)
        })
    })
    
    heatmap_layout["title"]["text"] = heatmapData[this_heatmap]["name"]
    
    Plotly.redraw("heatmap_space")
}

let heatmap_layout = {
    title: {
        text: ""
    },
    annotations: [],
    height: 600,
    width: 1200
}

let heatmap_plot = {
    z: [], 
    y: [],
    x: [],
    type: "heatmap",
    colorscale: "Viridis",
    annotations: []
}

let heatmap_buttons_to_remove = ["pan2d", "select2d", "lasso2d", "zoomIn2d", "zoomOut2d", "autoScale2d", "hoverClosestGl2d", "hoverClosestPie", "toggleHover", "toImage", "sendDataToCloud", "toggleSpikelines", "resetViewMapbox", "zoom"]

//On clicking a well in the heatmap, change the chroma_plot to show that particular well. 
function setHeatmapClickFns() {
    let myPlot = document.getElementById('heatmap_space')
    myPlot.on('plotly_click', function(clickdata){
        if(clickdata.points.length > 0 && clickdata.points[0].y !== undefined) {
            let well_id = clickdata.points[0].y + clickdata.points[0].x
            this_sample_no = samples.findIndex(x => x["well_readable"] == well_id)
            this_sample_name = samples[this_sample_no]["filename"]
            changeSample()

            //Stop showing the logs as we've switched to a "by-well" format. 
            $(' .logs ').addClass("d-none")

            //Set the name of the sample currently shown. 
            $(' #curr_sample_name ').html('Well: ' + samples[this_sample_no]["well_readable"] + '.')
        }
    });
}