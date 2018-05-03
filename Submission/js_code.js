// *************************************
// Zonotopes
// *************************************

// Functions for dragging interactivity
function dragstarted(d) {d3.select(this).raise().classed("active", true);}

function dragended(d) {d3.select(this).classed("active", false);}

function draggedA0Zono2D(d, i) {
    if (d[0] > 1) {
        A[i][0] = d[0] = 1;
    } else if (d[0] < -1) {
        A[i][0] = d[0] = -1;
    } else {
        A[i][0] = yAVecsZono2D.invert(d3.event.y);
    }
    
    A = MGS(A, 0);
    updateZonotope(A, m);
}

function draggedA1Zono2D(d, i) {
    if (d[1] > 1) {
        A[i][1] = d[1] = 1;
    } else if (d[1] < -1) {
        A[i][1] = d[1] = -1;
    } else {
        A[i][1] = yAVecsZono2D.invert(d3.event.y);
    }
    
    A = MGS(A, 1);
    updateZonotope(A, m);
}

// Function to compute the corners of a [-1,1] hypercube in m dimensions
function compute_hypercorners(mDim){
    var hc = new Array(2**mDim);
    for (i = 0; i < 2**mDim; i++) {
    hc[i] = new Array(mDim);
        for (j = 0; j < mDim; j++) {
            if (~~(i / (2**j) % 2) == 0) {
                hc[i][j] = -1;
            } else {
                hc[i][j] = 1;
            }
        }
    }
    return hc;
}

// Modified Gram-Schmidt function
function MGS(mA, index) {
    var mDim = mA.length;
    var nDim = mA[0].length;

    // transpose mA for computational convenience
    mA = numeric.transpose(mA);

    var indexList = d3.range(0, nDim, 1);
    
    var mA = mA.slice(indexList.indexOf(index), nDim)
               .concat(mA.slice(0, indexList.indexOf(index)));
    var indexList = indexList.slice(indexList.indexOf(index), nDim)
                             .concat(indexList.slice(0, indexList.indexOf(index)));

    // normalize first column of mA
   	var matNorm = numeric.norm2(mA[0]);
    for (i = 0; i < mDim; i++) {
    	mA[0][i] = mA[0][i]/matNorm;
    }

    // orthogonalize columns of mA
    for (i = 1; i < nDim; i++) {

    	for (j = 0; j < i; j++) {
    		var projUV = numeric.sum(numeric.mul(mA[i], mA[j]))/numeric.norm2Squared(mA[j]);
    		mA[i] = numeric.sub(mA[i], numeric.dot(mA[j], projUV));
    	}

		// normalize ith column of mA
	   	matNorm = numeric.norm2(mA[i]);
	    for (j = 0; j < mDim; j++) {
	    	mA[i][j] = mA[i][j]/matNorm;
	    }
    }
    
    var mA = mA.slice(indexList.indexOf(0), nDim)
               .concat(mA.slice(0, indexList.indexOf(0)));
    var indexList = indexList.slice(indexList.indexOf(0), nDim)
                             .concat(indexList.slice(0, indexList.indexOf(0)));

    return numeric.transpose(mA);
}

// Function for computing 2d convex hull of points
// Taken from: https://en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
function convexHull(points) {
    function cross(a, b, o) {
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    }  

    points.sort(function(a, b) {
        return a[0] == b[0] ? a[1] - b[1] : a[0] - b[0];
    });

    var lower = [];
    for (var i = 0; i < points.length; i++) {
        while (lower.length >= 2 && cross(lower[lower.length - 2], lower[lower.length - 1], points[i]) <= 0) {
            lower.pop();
        }
        lower.push(points[i]);
    }

    var upper = [];
    for (var i = points.length - 1; i >= 0; i--) {
        while (upper.length >= 2 && cross(upper[upper.length - 2], upper[upper.length - 1], points[i]) <= 0) {
            upper.pop();
        }
        upper.push(points[i]);
    }

    upper.pop();
    lower.pop();
    return lower.concat(upper);
}

// Function to compute the corners of the 2d zonotope
function compute_2d_zonotope(hc, mA, tf) {
    var mDim = hc[0].length;

    var zc = numeric.dot(hc, mA);

    var outerPoints = convexHull(zc);
    var zonoPoints = [];
    for (i = 0; i < outerPoints.length; i++) {
    	zonoPoints[i] = {data: outerPoints[i], class: 'OuterZonoPoint'};
    }

    if (tf) {
		for (i = 0; i < zc.length; i++) {
			if (outerPoints.indexOf(zc[i]) == -1) {
				zonoPoints.push({data: zc[i], class: 'InnerZonoPoint'});
			}
		}
	}
    return zonoPoints.reverse();
}

// Build svg for 2D zonotope plot
var marginZono2D = {top: 10, right:10, bottom: 50, left: 50},
    widthZono2D  = 500 - marginZono2D.left - marginZono2D.right,
    heightZono2D = 500 - marginZono2D.top - marginZono2D.bottom;

var svgZono2D = d3.select("#ZonotopePlots").append("svg")
                  .attr("width", widthZono2D + marginZono2D.left + marginZono2D.right)
                  .attr("height", heightZono2D + marginZono2D.top + marginZono2D.bottom)
                  .append("g")
                  .attr("transform", "translate(" + marginZono2D.left + "," + marginZono2D.top + ")");

// Build svg for plots of columns of A
var marginAVecsZono2D = {top: 30, right:10, bottom: 30, left: 30},
    widthAVecsZono2D  = 600 - marginAVecsZono2D.left - marginAVecsZono2D.right,
    heightAVecsZono2D = 300 - marginAVecsZono2D.top - marginAVecsZono2D.bottom;

var svgA0Zono2D = d3.select("#ZonotopeAPlots").append("svg")
                    .attr("width", widthAVecsZono2D + marginAVecsZono2D.left + marginAVecsZono2D.right)
                    .attr("height", heightAVecsZono2D + marginAVecsZono2D.top + marginAVecsZono2D.bottom)
                    .append("g")
                    .attr("transform", "translate(" + marginAVecsZono2D.left + "," + marginAVecsZono2D.top + ")");

d3.select("#ZonotopeAPlots").append("svg")
  .attr("width", 50)
  .attr("height", heightAVecsZono2D + marginAVecsZono2D.top + marginAVecsZono2D.bottom);

var svgA1Zono2D = d3.select("#ZonotopeAPlots").append("svg")
                    .attr("width", widthAVecsZono2D + marginAVecsZono2D.left + marginAVecsZono2D.right)
                    .attr("height", heightAVecsZono2D + marginAVecsZono2D.top + marginAVecsZono2D.bottom)
                    .append("g")
                    .attr("transform", "translate(" + marginAVecsZono2D.left + "," + marginAVecsZono2D.top + ")");

function mouseover2DZonotope(d, index, dataset) {
    d3.select(dataset[index])
      .transition()
      .duration(200)
      .attr("r", 10);
}

function mouseout2DZonotope(d, index, dataset) {
    d3.select(dataset[index])
      .transition()
      .duration(200)
      .attr("r", 5);
}

// Funciton to obtain new random rotation
function randomZonotope(mDim) {
	A = MGS(numeric.sub(numeric.random([mDim, 2]), 0.5), 0);
	
	updateZonotope(A, mDim);
}

// Function to update zonotope image
function updateZonotope(mA, mDim){

	var tfShowInnerPoints = [];
	d3.selectAll("#Zono2DCheckBox").each(function(){
          if (d3.select(this).property("checked")) {tfShowInnerPoints = true;}
          else {tfShowInnerPoints = false;}
        });

	// Compute corners of zonotope
	var zonotopeCorners = compute_2d_zonotope(compute_hypercorners(mDim), mA, tfShowInnerPoints);

	// Scale the range of the data
	xZono2D.domain([-Math.sqrt(mDim), Math.sqrt(mDim)]);
	yZono2D.domain([-Math.sqrt(mDim), Math.sqrt(mDim)]);

	// Plot zonotope shape
	svgZono2D.selectAll("polygon")
	         .remove();
	svgZono2D.selectAll("polygon")
	         .data([zonotopeCorners])
	         .enter()
	         .append("polygon")
	         .attr("class", "Zonotope")
	         .attr("points", function(d) { 
	         	  return d.map(function(d){
	         		  if (d.class == "OuterZonoPoint") {
	         		  	  return [xZono2D(d.data[0]),yZono2D(d.data[1])].join(",");
	         		  }
	         	  }).join(" ");
	    	  });

	// Add the scatterplot points
	svgZono2D.selectAll("circle")
	         .remove();
	svgZono2D.selectAll("circle")
	         .data(zonotopeCorners)
	         .enter()
	         .append("circle")
	         .attr("class", function(d) {return d.class;})
	         .attr("r", 5)
	         .attr("cx", function(d) {return xZono2D(d.data[0]);})
	         .attr("cy", function(d) {return yZono2D(d.data[1]);})
           .on("mouseover", mouseover2DZonotope)
           .on("mouseout", mouseout2DZonotope);
    
  // Update the axes
  svgZono2D.select("#xAxisZono2D").call(xAxisZono2D);
  svgZono2D.select("#yAxisZono2D").call(yAxisZono2D);

	// Scale the range of the data
	xAVecsZono2D.domain([0, mDim+1]);
	yAVecsZono2D.domain([-1, 1]);

	// Add the points for entries of A
	svgA0Zono2D.selectAll("circle")
	           .remove();
	svgA0Zono2D.selectAll("circle")
	           .data(mA)
	           .enter()
	           .append("circle")
	           .attr("class", "Acomponents")
	           .attr("r", 10)
	           .attr("cx", function(d, i) {return xAVecsZono2D(i+1);})
	           .attr("cy", function(d) {return yAVecsZono2D(d[0]);})
             .call(d3.drag()
                  .on("start", dragstarted)
                  .on("drag", draggedA0Zono2D)
                  .on("end", dragended));

  // Update the axes
  svgA0Zono2D.select("#xAxisAVecsZono2D").call(xAxisAVecsZono2D);
  svgA0Zono2D.select("#yAxisAVecsZono2D").call(yAxisAVecsZono2D);

	svgA1Zono2D.selectAll("circle")
	           .remove();
	svgA1Zono2D.selectAll("circle")
	           .data(mA)
	           .enter()
	           .append("circle")
	           .attr("class", "Acomponents")
	           .attr("r", 10)
	           .attr("cx", function(d, i) {return xAVecsZono2D(i+1);})
	           .attr("cy", function(d) {return yAVecsZono2D(d[1]);})
             .call(d3.drag()
                  .on("start", dragstarted)
                  .on("drag", draggedA1Zono2D)
                  .on("end", dragended));

    // Update the axes
    svgA1Zono2D.select("#xAxisAVecsZono2D").call(xAxisAVecsZono2D);
    svgA1Zono2D.select("#yAxisAVecsZono2D").call(yAxisAVecsZono2D);
}

// Initialize m and A
var m = 3;
var A = MGS(numeric.sub(numeric.random([m, 2]), 0.5), 0);

var xZono2D = d3.scaleLinear().domain([-Math.sqrt(m), Math.sqrt(m)]).range([0, widthZono2D]),
    yZono2D = d3.scaleLinear().domain([-Math.sqrt(m), Math.sqrt(m)]).range([heightZono2D, 0]),
    xAxisZono2D = d3.axisBottom(xZono2D).tickValues([-20, 20]),
    yAxisZono2D = d3.axisLeft(yZono2D).tickValues([-20, 20]);


	// Scale the range of the data
var xAVecsZono2D = d3.scaleLinear().domain([0, m+1]).range([0, widthAVecsZono2D]),
    yAVecsZono2D = d3.scaleLinear().domain([-1, 1]).range([heightAVecsZono2D, 0]),
    xAxisAVecsZono2D = d3.axisBottom(xAVecsZono2D).tickValues([-20, 20]), //.tickValues(d3.range(1, m+1, 1)),
    yAxisAVecsZono2D = d3.axisLeft(yAVecsZono2D);

// Add the X & Y Axes to zonotope plot
svgZono2D.append("g")
         .attr("transform", "translate(0, " + heightZono2D + ")")
         .call(xAxisZono2D);
svgZono2D.append("g").call(yAxisZono2D);


// Add the X & Y Axes to a1 & a2 vector plots
svgA0Zono2D.append("g")
           .attr("transform", "translate(0, " + heightAVecsZono2D / 2 + ")")
           .call(xAxisAVecsZono2D);
svgA0Zono2D.append("g").call(yAxisAVecsZono2D);
svgA1Zono2D.append("g")
           .attr("transform", "translate(0, " + heightAVecsZono2D / 2 + ")")
           .call(xAxisAVecsZono2D);
svgA1Zono2D.append("g").call(yAxisAVecsZono2D);

// Set callback functions
d3.selectAll("#Zono2DCheckBox")
  .on("change", function(){updateZonotope(A, m);});
d3.selectAll("#Zono2DRefreshButton")
  .on("click", function(){randomZonotope(m);});
d3.selectAll("#Zono2DSelectDimension")
  .on("change", function(a, b, c) {m = c[0].selectedIndex + 3; randomZonotope(m);});

// Place axes labels
svgZono2D.append("text")
         .attr("class", "label")
         .text("a1^T x")
         .attr("x", widthZono2D/2)
         .attr("y", heightZono2D + 40);
svgZono2D.append("text")
         .attr("class", "label")
         .text("a2^T x")
         .attr("x", heightZono2D/2)
         .attr("y", 35)
         .attr("transform", "rotate(90)")
         .style("text-anchor", "start");
svgA0Zono2D.append("text")
           .attr("class", "label")
           .text("a1")
           .attr("x", (widthAVecsZono2D / 2))
           .attr("y", 0 - (marginAVecsZono2D.top / 2));
svgA1Zono2D.append("text")
           .attr("class", "label")
           .text("a2")
           .attr("x", (widthAVecsZono2D / 2))
           .attr("y", 0 - (marginAVecsZono2D.top / 2));

updateZonotope(A, m);

// *************************************
// 3D Zonotopes?!?!?!
// *************************************

/*
var svgZono3D = d3.select("#ZonotopePlots").append("svg")
                  .attr("width", 100)
                  .attr("height", 100)
                  .append("g")
                  .attr("transform", "translate(10, 120)");
*/

// *************************************
// Active Subspaces
// *************************************

// Evaluate the piston function
// Adapted from: https://www.sfu.ca/~ssurjano/piston.html
function pistonFunction(xx) {
    var MM = xx.length;

    xx = numeric.transpose(xx);

    // Scale inputs to physically meaningful ranges
    xx[0] = rescaleInput(xx[0], 30, 60); // M
    xx[1] = rescaleInput(xx[1], 0.005, 0.02); // S
    xx[2] = rescaleInput(xx[2], 0.002, 0.01); // V0
    xx[3] = rescaleInput(xx[3], 1000, 5000); // k
    xx[4] = rescaleInput(xx[4], 90000, 111000); // P0
    xx[5] = rescaleInput(xx[5], 290, 296); // Ta
    xx[6] = rescaleInput(xx[6], 340, 360); // T0

    var Aterm1 = numeric.mul(xx[4], xx[1]),
        Aterm2 = numeric.mul(19.62, xx[0]),
        Aterm3 = numeric.div(numeric.mul(xx[3], xx[2]), xx[1]);
    var Aterm  = numeric.sub(numeric.add(Aterm1, Aterm2), Aterm3);

    var Vsubterm1 = numeric.pow(Aterm, 2),
        Vsubterm2 = numeric.div(numeric.mul(xx[4], xx[2]), xx[6]),
        Vsubterm3 = numeric.mul(4, numeric.mul(xx[3], numeric.mul(xx[5])));

    var Vterm1 = numeric.div(xx[1], numeric.mul(2, xx[3])),
        Vterm2 = numeric.add(Vsubterm1, numeric.mul(Vsubterm2, Vsubterm3));
    var Vterm  = numeric.mul(Vterm1, numeric.sub(numeric.sqrt(Vterm2), Aterm));

    var Cterm1 = numeric.pow(xx[1], 2),
        Cterm2 = numeric.div(xx[5], numeric.pow(Vterm, 2));
    var Cterm  = numeric.mul(Cterm1, numeric.mul(Vsubterm2, Cterm2));

    var ff = numeric.mul(2, numeric.sqrt(numeric.div(xx[0], numeric.add(xx[3], Cterm))));

    return numeric.transpose([ff]);
}

function rescaleInput(xx, lb, ub) {
  return numeric.add(numeric.mul((ub-lb)/2, numeric.add(xx, 1)), lb);
}

function generateInputs(MM) {
  var xx = new Array(MM);
  for (i = 0; i < MM; i++) {
    var xxTemp = numeric.sub(numeric.mul(2, numeric.random([m2])), 1);
    
    if (numeric.random([1]) < Math.min(...numeric.abs(xxTemp))/100) {
      xx[i] = xxTemp;
    } else {
      i--;
    }
  }
  return xx;
}

function zipData(xx, ff) {
  var MM = xx.length,
      data = new Array(MM);
  for (i = 0; i < MM; i++) {
    data[i] = {x0: xx[i][0], x1: xx[i][1], f: ff[i]};
  }
  return data;
}

function resetAS(mW) {
  W = mW;
  updateActiveSubspaces(X, f, W);
}

function draggedA0AS(d, i) {
  if (d[0] > 1) {
      W[i][0] = d[0] = 1;
  } else if (d[0] < -1) {
      W[i][0] = d[0] = -1;
  } else {
      W[i][0] = yAVecsAS.invert(d3.event.y);
  }
    
  W = MGS(W, 0);
  updateActiveSubspaces(X, f, W);
}

function draggedA1AS(d, i) {
  if (d[1] > 1) {
      W[i][1] = d[1] = 1;
  } else if (d[1] < -1) {
      W[i][1] = d[1] = -1;
  } else {
      W[i][1] = yAVecsAS.invert(d3.event.y);
  }
    
  W = MGS(W, 1);
  updateActiveSubspaces(X, f, W);
}

// Build svg for active subspace plot
var marginAS = {top: 10, right:10, bottom: 50, left: 50},
    widthAS  = 500 - marginAS.left - marginAS.right,
    heightAS = 500 - marginAS.top - marginAS.bottom;

var svgAS1D = d3.select("#activeSubspaces").append("svg")
                .attr("width", widthAS + marginAS.left + marginAS.right)
                .attr("height", heightAS + marginAS.top + marginAS.bottom)
                .append("g")
                .attr("transform", "translate(" + marginAS.left + "," + marginAS.top + ")");

d3.select("#activeSubspaces").append("svg")
  .attr("width", 50)
  .attr("height", heightAS + marginAS.top + marginAS.bottom);

var svgAS2D = d3.select("#activeSubspaces").append("svg")
                .attr("width", widthAS + marginAS.left + marginAS.right)
                .attr("height", heightAS + marginAS.top + marginAS.bottom)
                .append("g")
                .attr("transform", "translate(" + marginAS.left + "," + marginAS.top + ")");

// Build svg for plots of columns of A
var marginAVecsAS = {top: 30, right:10, bottom: 30, left: 40},
    widthAVecsAS  = 600 - marginAVecsAS.left - marginAVecsAS.right,
    heightAVecsAS = 300 - marginAVecsAS.top - marginAVecsAS.bottom;

var svgA0AS = d3.select("#activeSubspacesAPlots").append("svg")
                .attr("width", widthAVecsAS + marginAVecsAS.left + marginAVecsAS.right)
                .attr("height", heightAVecsAS + marginAVecsAS.top + marginAVecsAS.bottom)
                .append("g")
                .attr("transform", "translate(" + marginAVecsAS.left + "," + marginAVecsAS.top + ")");

d3.select("#activeSubspacesAPlots").append("svg")
  .attr("width", 50)
  .attr("height", heightAVecsAS + marginAVecsAS.top + marginAVecsAS.bottom);

var svgA1AS = d3.select("#activeSubspacesAPlots").append("svg")
                .attr("width", widthAVecsAS + marginAVecsAS.left + marginAVecsAS.right)
                .attr("height", heightAVecsAS + marginAVecsAS.top + marginAVecsAS.bottom)
                .append("g")
                .attr("transform", "translate(" + marginAVecsAS.left + "," + marginAVecsAS.top + ")");


// Function to update the active subspaces plots
function updateActiveSubspaces(xx, ff, mW){

  var data = zipData(numeric.dot(xx, mW), ff);

  // Add the scatterplot points
  svgAS1D.selectAll("circle")
         .remove();
  svgAS1D.selectAll("circle")
         .data(data)
         .enter()
         .append("circle")
         .attr("class", "OneDShadowPlot")
         .attr("r", 5)
         .attr("cx", function(d) {return xAS1D(d.x0);})
         .attr("cy", function(d) {return yAS1D(d.f);})
         .on("mouseover", mouseover2DZonotope)
         .on("mouseout", mouseout2DZonotope);

  // Add the scatterplot points
  svgAS2D.selectAll("circle")
         .remove();
  svgAS2D.selectAll("circle")
         .data(data)
         .enter()
         .append("circle")
         .attr("class", "TwoDShadowPlot")
         .attr("r", 5)
         .attr("cx", function(d) {return xAS2D(d.x0);})
         .attr("cy", function(d) {return yAS2D(d.x1);})
         .attr("fill", function(d) {return cAS2D(d.f);})
         .on("mouseover", mouseover2DZonotope)
         .on("mouseout", mouseout2DZonotope);
  
  // Add the points for entries of A
  svgA0AS.selectAll("circle")
         .remove();
  svgA0AS.selectAll("circle")
         .data(mW)
         .enter()
         .append("circle")
         .attr("class", "Acomponents")
         .attr("r", 10)
         .attr("cx", function(d, i) {return xAVecsAS(i+1);})
         .attr("cy", function(d) {return yAVecsAS(d[0]);})
         .call(d3.drag()
              .on("start", dragstarted)
              .on("drag", draggedA0AS)
              .on("end", dragended));

  svgA1AS.selectAll("circle")
         .remove();
  svgA1AS.selectAll("circle")
         .data(mW)
         .enter()
         .append("circle")
         .attr("class", "Acomponents")
         .attr("r", 10)
         .attr("cx", function(d, i) {return xAVecsAS(i+1);})
         .attr("cy", function(d) {return yAVecsAS(d[1]);})
         .call(d3.drag()
              .on("start", dragstarted)
              .on("drag", draggedA1AS)
              .on("end", dragended));

}

// Set the scale for the 1D and 2D shadow plots
var xAS1D = d3.scaleLinear().domain([-Math.sqrt(7), Math.sqrt(7)]).range([0, widthAS]),
    yAS1D = d3.scaleLinear().domain([0, 0.4]).range([heightAS, 0]),
    xAxisAS1D = d3.axisBottom(xAS1D),
    yAxisAS1D = d3.axisLeft(yAS1D);

var xAS2D = d3.scaleLinear().domain([-Math.sqrt(7), Math.sqrt(7)]).range([0, widthAS]),
    yAS2D = d3.scaleLinear().domain([-Math.sqrt(7), Math.sqrt(7)]).range([heightAS, 0]),
    cAS2D = d3.scaleLinear().domain([0.05, 0.3]).range(["#CFEBFF", "#0000FF"]),
    //cAS2D = d3.scaleQuantile().domain([0.05, 0.3]).range(colorbrewer.YlGnBu[7]),
    xAxisAS2D = d3.axisBottom(xAS2D),
    yAxisAS2D = d3.axisLeft(yAS2D);
    //cAS2D = d3.scaleLinear().domain([0, 0.4]).range(["#87CEFF", "#0000FF"]),

// Set the scale for the vectors of A plots
var xAVecsAS = d3.scaleLinear().domain([0, 8]).range([0, widthAVecsAS]),
    yAVecsAS = d3.scaleLinear().domain([-1, 1]).range([heightAVecsAS, 0]),
    xAxisAVecsAS = d3.axisBottom(xAVecsAS).tickValues(d3.range(1, 8, 1)),
    yAxisAVecsAS = d3.axisLeft(yAVecsAS);

// Add the X & Y Axes to the 1D and 2D shadow plots
svgAS1D.append("g")
       .attr("transform", "translate(0, " + heightAS + ")")
       .call(xAxisAS1D);
svgAS1D.append("g").call(yAxisAS1D);
svgAS2D.append("g")
       .attr("transform", "translate(0, " + heightAS + ")")
       .call(xAxisAS2D);
svgAS2D.append("g").call(yAxisAS2D);


// Add the X & Y Axes to a1 & a2 vector plots
svgA0AS.append("g")
       .attr("transform", "translate(0, " + heightAVecsAS / 2 + ")")
       .call(xAxisAVecsAS);
svgA0AS.append("g").call(yAxisAVecsAS);
svgA1AS.append("g")
       .attr("transform", "translate(0, " + heightAVecsAS / 2 + ")")
       .call(xAxisAVecsAS);
svgA1AS.append("g").call(yAxisAVecsAS);


// Add axes labels to the 1D and 2D shadow plots
svgAS1D.append("text")
       .attr("class", "label")
       .text("a1^T x")
       .attr("x", widthAS/2)
       .attr("y", heightAS + 40);
svgAS1D.append("text")
       .attr("class", "label")
       .text("y")
       .attr("x", heightAS/2)
       .attr("y", 45)
       .attr("transform", "rotate(90)")
       .style("text-anchor", "start");
svgAS2D.append("text")
       .attr("class", "label")
       .text("a1^T x")
       .attr("x", widthAS/2)
       .attr("y", heightAS + 40);
svgAS2D.append("text")
       .attr("class", "label")
       .text("a2^T x")
       .attr("x", heightAS/2)
       .attr("y", 45)
       .attr("transform", "rotate(90)")
       .style("text-anchor", "start");

// Add axes labesl to a1 & a2 vector plots
svgA0AS.append("text")
       .attr("class", "label")
       .text("a1")
       .attr("x", (widthAVecsAS / 2))
       .attr("y", 0 - (marginAVecsAS.top / 2));
svgA1AS.append("text")
       .attr("class", "label")
       .text("a2")
       .attr("x", (widthAVecsAS / 2))
       .attr("y", 0 - (marginAVecsAS.top / 2));


d3.selectAll("#ASResetButton")
  .on("click", function(){resetAS(W_AS);});

var W_AS = [[ 0.160602933911370,   0.335264496144912],
            [-0.794311961306026,  -0.319519481384172],
            [ 0.575703067832273,  -0.645049603954946],
            [-0.104388344185824,  -0.607478171762296],
            [-0.030406736335403,  -0.012146886452726],
            [ 0.001496577237633,  -0.005212322904076],
            [-0.004177041284026,   0.014545884299044]];

var W = [[ 0.160602933911370,   0.335264496144912],
         [-0.794311961306026,  -0.319519481384172],
         [ 0.575703067832273,  -0.645049603954946],
         [-0.104388344185824,  -0.607478171762296],
         [-0.030406736335403,  -0.012146886452726],
         [ 0.001496577237633,  -0.005212322904076],
         [-0.004177041284026,   0.014545884299044]];

// Initialize m and A
var m2 = 7,
    M  = 5e2;
var X = generateInputs(M);
//var X = numeric.sub(numeric.mul(2, numeric.random([M, m2])), 1);

var f = pistonFunction(X);

updateActiveSubspaces(X, f, W);

// Get key directions from the piston model (use MATLAB)
// ---hardcode those directions and have a reset button set the AS
// ---also button with random directions
// Make 1d and 2d shadow plots as well as the 2 vectors of A defining the reduction
// Give user control to rotate the space same as before

