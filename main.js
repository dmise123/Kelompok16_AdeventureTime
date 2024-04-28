function generateMount(num_triangles, radius, height, topColor, bodyColor, scale) {
    var cone_vertices = [];
    var cone_colors = [];
    var angle_increment = (2 * Math.PI) / num_triangles;
  
  
    // Generate cone vertices and colors for the body
    for (var i = 0; i < num_triangles; i++) {
        var angle1 = i * angle_increment;
        var angle2 = (i + 1) * angle_increment;
  
  
        // Apex vertex
        cone_vertices.push(0, height, 0);
        cone_colors.push(...topColor); // Set top segment color
  
  
        // First base vertex
        cone_vertices.push(radius * Math.cos(angle1) * scale, 0, radius * Math.sin(angle1) * scale);
        cone_colors.push(...bodyColor); // Set body segment color
  
  
        // Second base vertex
        cone_vertices.push(radius * Math.cos(angle2) * scale, 0, radius * Math.sin(angle2) * scale);
        cone_colors.push(...bodyColor); // Set body segment color
    }
  
  
    // Define vertices for the bottom of the cone (circle)
    var bottom_center = [0, 0, 0];
    cone_vertices.push(...bottom_center); // Center vertex for bottom
    cone_colors.push(...bodyColor); // Set body segment color for bottom center
  
  
    // Generate indices for rendering the base of the cone
    for (var i = 0; i < num_triangles; i++) {
        var angle1 = i * angle_increment;
        var angle2 = (i + 1) * angle_increment;
  
  
        cone_vertices.push(radius * Math.cos(angle2) * scale, 0, radius * Math.sin(angle2) * scale); // Vertex on base
        cone_colors.push(...bodyColor); // Set body segment color for base vertex
  
  
        cone_vertices.push(radius * Math.cos(angle1) * scale, 0, radius * Math.sin(angle1) * scale); // Next vertex on base
        cone_colors.push(...bodyColor); // Set body segment color for base vertex
  
  
        // Connecting vertex between base and center
        cone_vertices.push(0, 0, 0); // Center vertex
        cone_colors.push(...bodyColor); // Set body segment color for center vertex
    }
  
  
    return {
        vertices: cone_vertices,
        colors: cone_colors
    };
  }
  function generateConnectedZigzagPath(numSegments, width, height, numZigs, scale ) {
    const vertices = [];
    const indices = [];
    const colors = [];
  
  
    var decrement = 1;
    var increment = 0;
  
  
    const segmentHeight = height / numSegments;
    const zigWidth = width / (numZigs * 2);
  
  
    // Generate vertices and colors for zigzag paths
    for (let i = 0; i <= numSegments; i++) {
        const y = i * segmentHeight;
        const isEvenRow = i % 2 === 0;
  
  
        // Generate zigzag pattern within the width
        for (let j = 0; j <= numZigs; j++) {
          
            const x = isEvenRow ? j * zigWidth : width - j * zigWidth;
            const z = Math.random() * segmentHeight; // Random elevation for each zigzag vertex
  
  
            vertices.push((x - width / 2) * decrement * scale, y  * scale , z * scale * increment); // Push vertices centered around the origin
            decrement -= 0.01;
            increment += 0.04;
            // Assign color to the vertex
            const brightness = 0.8 + 0.2 * Math.random(); // Random brightness between 0.8 and 1.0
            colors.push(0.1, 207/255 * Math.random(), 249/255); // Grayish color
        }
    }
  
  
    // Generate indices for rendering triangles
    for (let i = 0; i < numSegments; i++) {
        const startIndex = i * (numZigs + 1);
        for (let j = 0; j < numZigs; j++) {
            // Define indices for two triangles forming a quad
            const index1 = startIndex + j;
            const index2 = startIndex + j + 1;
            const index3 = startIndex + j + numZigs + 1;
            const index4 = startIndex + j + numZigs + 2;
  
  
            // Push indices for the first triangle
            indices.push(index1, index2, index3);
            // Push indices for the second triangle
            indices.push(index2, index4, index3);
  
  
            // Connect to the corresponding vertex of the zigzag above
            if (i < numSegments - 1) {
                const nextStartIndex = (i + 1) * (numZigs + 1);
                const topIndex1 = nextStartIndex + j;
                const topIndex2 = nextStartIndex + j + 1;
  
  
                // Connect the current vertex to the corresponding top vertex
                indices.push(index2, topIndex2, index4);
                indices.push(index2, topIndex1, topIndex2);
            }
        }
    }
  
  
    return {
        vertices: vertices,
        indices: indices,
        colors: colors // Include colors array in the returned object
    };
  }
  
  
  function generateHalfEllipsoid(radius, segments, sphereColor, baseColor, topColor) {
    const vertices = [];
    const indices = [];
    const colors = [];
  
  
    // Generate vertices for the half sphere
    for (let i = 0; i <= segments / 2; i++) {
        const lat = Math.PI * i / segments;
        const sinLat = Math.sin(lat);
        const cosLat = Math.cos(lat);
  
  
        for (let j = 0; j <= segments; j++) {
            const lng = 2 * Math.PI * j / segments;
            const x = cosLat * Math.cos(lng);
            const y = cosLat * Math.sin(lng);
            const z = sinLat;
  
  
            vertices.push(radius * x, radius * y, radius * z);
  
  
            // Assign sphere color or top color based on vertex position
            if (i < segments / 4) {
                colors.push(...sphereColor); // Top part of the half sphere
            } else {
                colors.push(...topColor); // Bottom part of the half sphere
            }
        }
    }
  
  
    // Add vertices for the base circle
    for (let j = 0; j <= segments; j++) {
        const lng = 2 * Math.PI * j / segments;
        const x = radius * Math.cos(lng);
        const y = radius * Math.sin(lng);
        const z = 0; // Base at the lowest point of the half sphere
  
  
        vertices.push(x, y, z);
        colors.push(...baseColor);
    }
  
  
    // Indices for the spherical part
    for (let i = 0; i < segments / 2; i++) {
        for (let j = 0; j < segments; j++) {
            const first = i * (segments + 1) + j;
            const second = first + segments + 1;
  
  
            indices.push(first, second, first + 1);
            indices.push(second, second + 1, first + 1);
        }
    }
  
  
    // Indices for the base circle
    const baseCenterIndex = vertices.length / 3; // Index of the center vertex for the base
    const firstBaseVertIndex = baseCenterIndex - segments - 1; // Index of the first vertex on the base perimeter
    for (let j = 0; j < segments; j++) {
        const current = firstBaseVertIndex + j;
        const next = current + 1;
  
  
        indices.push(baseCenterIndex, current, next);
    }
  
  
    return { vertices, indices, colors };
  }
  function generateElipticPara(radius, segments, scale, r, g, b, condition) {
    const vertices = [];
    const indices = [];
    const colors = [];
  
  
    const segmentColors = [
        [0, 1, 0.0],
        [0, 0.8, 0.0], // Gray
        [0, 0.6, 0.0],
        [0, 0.4, 0], // Yellow
    ];
  
  
  
  
    for (let i = 0; i <= segments; i++) {
        const lat = Math.PI / 2 * i / segments;
        const xy = radius * Math.sqrt(Math.pow(Math.sin(lat), 2));
        var z = Math.pow(lat, 2) * xy * scale;
  
  
        for (let j = 0; j <= segments; j++) {
            const lng = 2 * Math.PI * j / segments;
            var x = xy * Math.cos(lng) * scale;
            var y = xy * Math.sin(lng) * scale;
  
  
            if (condition == true) {
                vertices.push(x, y, z - 1.5 * Math.random());
                let colorIndex;
                if (i <= 82) colorIndex = 0; // Blue for segments 1-5
                else if (i > 82 && i <= 92) colorIndex = 1; // Gray for segments 6-7
                else if (i > 92 && i <= 97) colorIndex = 2; // Blue for segment 8
  
  
                else colorIndex = 3; // Dark Gray for segment 10
  
  
                colors.push(...segmentColors[colorIndex]);
            } else {
                vertices.push(x, y, z - 1.5);
                colors.push(r, g, b);
            }
  
  
  
  
            // Assign color based on segment index
  
  
  
  
  
  
        }
    }
  
  
    for (let i = 0; i < segments; i++) {
        for (let j = 0; j < segments; j++) {
            const first = i * (segments + 1) + j;
            const second = first + segments + 1;
  
  
            indices.push(first, second, first + 1);
            indices.push(second, second + 1, first + 1);
        }
    }
  
  
    return { vertices, indices, colors };
  }
  function halfbezierTorus(mainRadius, tubeRadius, numMainPoints, numTubePoints, scale, color) {
    var mainControlPoints = [];
    var tubeControlPoints = [];
    var mainAngleIncrement = Math.PI / (numMainPoints - 1); // Only half circle for main circle
    var tubeAngleIncrement = (2 * Math.PI) / (numTubePoints - 1);
  
  
    // Create control points along the main half circle of the torus
    for (var i = 0; i < numMainPoints; i++) {
        var mainAngle = i * mainAngleIncrement;
        var mainX = mainRadius * Math.cos(mainAngle);
        var mainY = mainRadius * Math.sin(mainAngle);
        var mainZ = 0; // Torus is centered at the origin for simplicity
  
  
        mainControlPoints.push({ x: mainX, y: mainY, z: mainZ });
    }
  
  
    // Create control points along the full circle for the tube
    for (var j = 0; j < numTubePoints; j++) {
        var tubeAngle = j * tubeAngleIncrement;
        var tubeX = tubeRadius * Math.cos(tubeAngle);
        var tubeY = tubeRadius * Math.sin(tubeAngle);
        var tubeZ = 0; // Torus tube is centered around the tube circle center
  
  
        tubeControlPoints.push({ x: tubeX, y: tubeY, z: tubeZ });
    }
  
  
    // Interpolate and generate bezier curve points for the main and tube
    var mainCurvePoints = bezierCurve(mainControlPoints, numMainPoints);
    var tubeCurvePoints = bezierCurve(tubeControlPoints, numTubePoints);
  
  
    // Combine main and tube curves to form the half-torus
    var torusPoints = [];
    var torusColors = [];
    var torusIndices = [];
  
  
    // Generate vertices and assign colors
    for (var k = 0; k < numMainPoints; k++) {
        for (var l = 0; l < numTubePoints; l++) {
            var torusPoint = {
                x: (mainCurvePoints[k].x + tubeCurvePoints[l].x) * scale,
                y: (mainCurvePoints[k].y + tubeCurvePoints[l].y) * scale,
                z: (mainCurvePoints[k].z + tubeCurvePoints[l].z) * scale
            };
  
  
            torusPoints.push(torusPoint.x, torusPoint.y, torusPoint.z);
            torusColors.push(...color); // Yellow color for demonstration
        }
    }
  
  
    // Generate indices to form triangles
    for (var m = 0; m < numMainPoints - 1; m++) {
        for (var n = 0; n < numTubePoints - 1; n++) {
            var currentIndex = m * numTubePoints + n;
            var nextIndex = currentIndex + numTubePoints;
  
  
            torusIndices.push(currentIndex, nextIndex, currentIndex + 1);
            torusIndices.push(nextIndex, nextIndex + 1, currentIndex + 1);
        }
    }
  
  
    return { vertices: torusPoints, indices: torusIndices, colors: torusColors };
  }
  function bezierTorus(mainRadius, tubeRadius, numMainPoints, numTubePoints, scale, color) {
    var mainControlPoints = [];
    var tubeControlPoints = [];
    var mainAngleIncrement = (2 * Math.PI) / (numMainPoints - 1);
    var tubeAngleIncrement = (2 * Math.PI) / (numTubePoints - 1);
  
  
    // Create control points along the main circle of the torus
    for (var i = 0; i < numMainPoints; i++) {
        var mainAngle = i * mainAngleIncrement;
        var mainX = mainRadius * Math.cos(mainAngle);
        var mainY = mainRadius * Math.sin(mainAngle);
        var mainZ = 0; // Torus is centered at the origin for simplicity
  
  
        mainControlPoints.push({ x: mainX, y: mainY, z: mainZ });
    }
  
  
    // Create control points along the tube of the torus
    for (var j = 0; j < numTubePoints; j++) {
        var tubeAngle = j * tubeAngleIncrement;
        var tubeX = (mainRadius + tubeRadius * Math.cos(tubeAngle)) * Math.cos(mainAngle);
        var tubeY = (mainRadius + tubeRadius * Math.cos(tubeAngle)) * Math.sin(mainAngle);
        var tubeZ = tubeRadius * Math.sin(tubeAngle);
  
  
        tubeControlPoints.push({ x: tubeX, y: tubeY, z: tubeZ });
    }
  
  
    // Interpolate and generate bezier curve points for the torus
    var mainCurvePoints = bezierCurve(mainControlPoints, numMainPoints);
    var tubeCurvePoints = bezierCurve(tubeControlPoints, numTubePoints);
  
  
    // Combine main and tube curves to form the torus
    var torusPoints = [];
    var torusColors = [];
    var torusIndices = [];
  
  
    // Generate vertices and assign colors
    for (var k = 0; k < numMainPoints; k++) {
        for (var l = 0; l < numTubePoints; l++) {
            var torusPoint = {
                x: (mainCurvePoints[k].x + tubeCurvePoints[l].x) * scale,
                y: (mainCurvePoints[k].y + tubeCurvePoints[l].y) * scale,
                z: (mainCurvePoints[k].z + tubeCurvePoints[l].z) * scale
            };
  
  
            torusPoints.push(torusPoint.x, torusPoint.y, torusPoint.z);
  
  
            // Assign color based on segment index or other criteria
            // Example: Assigning a uniform color for demonstration
            torusColors.push(...color); // Yellow color for demonstration
        }
    }
  
  
    // Generate indices to form triangles
    for (var m = 0; m < numMainPoints - 1; m++) {
        for (var n = 0; n < numTubePoints - 1; n++) {
            var currentIndex = m * numTubePoints + n;
            var nextIndex = currentIndex + numTubePoints;
  
  
            torusIndices.push(currentIndex, nextIndex, currentIndex + 1);
            torusIndices.push(nextIndex, nextIndex + 1, currentIndex + 1);
        }
    }
  
  
    return { vertices: torusPoints, indices: torusIndices, colors: torusColors };
  }
  function bezierCurve(controlPoints, numPoints) {
    var curvePoints = [];
    var n = controlPoints.length - 1;
  
  
    for (var i = 0; i < numPoints; i++) {
        var t = i / (numPoints - 1);
        var point = { x: 0, y: 0, z: 0 };
  
  
        for (var j = 0; j <= n; j++) {
            var b = binomialCoefficient(n, j) * Math.pow(1 - t, n - j) * Math.pow(t, j);
            point.x += controlPoints[j].x * b;
            point.y += controlPoints[j].y * b;
            point.z += controlPoints[j].z * b;
        }
  
  
        curvePoints.push(point);
    }
  
  
    return curvePoints;
  }
  function binomialCoefficient(n, k) {
    if (k === 0 || k === n) {
        return 1;
    }
  
  
    var coeff = 1;
    for (var x = n - k + 1; x <= n; x++) {
        coeff *= x;
    }
    for (x = 1; x <= k; x++) {
        coeff /= x;
    }
  
  
    return coeff;
  }
  function generateEllipsoid(x, y, z, a, b, c, segments, color) {
    var vertices = [];
    var indices = [];
    var colors = [];
  
  
    for (var j = 0; j <= segments; j++) {
        var theta = j * Math.PI / segments;
        var sinTheta = Math.sin(theta);
        var cosTheta = Math.cos(theta);
  
  
        for (var i = 0; i <= segments; i++) {
            var phi = i * 2 * Math.PI / segments;
            var sinPhi = Math.sin(phi);
            var cosPhi = Math.cos(phi);
  
  
            var vx = x + a * cosPhi * sinTheta;
            var vy = y + b * cosTheta;
            var vz = z + c * sinPhi * sinTheta;
  
  
            vertices.push(vx, vy - 0.7, vz);
  
  
            // Assign color to each vertex
            colors.push(...color);
        }
    }
  
  
    for (var j = 0; j < segments; j++) {
        var rowStart = j * (segments + 1);
        var nextRowStart = (j + 1) * (segments + 1);
  
  
        for (var i = 0; i < segments; i++) {
            indices.push(rowStart + i, nextRowStart + i);
        }
  
  
        indices.push(rowStart, nextRowStart);
    }
  
  
    return {
        vertices: vertices,
        indices: indices,
        colors: colors
    };
  }
  function createSphere(radius, segments) {
    const vertices = [];
    const indices = [];
    const colors = [];
  
  
    // Define colors for each segment
    const segmentColors = [
        [0.5, 0.5, 1.0],
        [0.8, 0.8, 0.8], // Gray
        [0.5, 0.5, 1.0],
        [1, 1, 0], // Yellow
        [0.4, 0.4, 0.4] // Dark Gray
    ];
  
  
    for (let i = 0; i <= segments; i++) {
        const lat = Math.PI * i / segments;
        const sinLat = Math.sin(lat);
        const cosLat = Math.cos(lat);
  
  
        for (let j = 0; j <= segments; j++) {
            const lng = 2 * Math.PI * j / segments;
            const sinLng = Math.sin(lng);
            const cosLng = Math.cos(lng);
  
  
            const x = cosLng * sinLat;
            const y = cosLat;
            const z = sinLng * sinLat;
  
  
            vertices.push(radius * x, radius * y, radius * z);
  
  
            // Assign color based on segment index
            let colorIndex;
            if (i <= 82) colorIndex = 0; // Blue for segments 1-5
            else if (i > 82 && i <= 92) colorIndex = 1; // Gray for segments 6-7
            else if (i > 92 && i <= 97) colorIndex = 0; // Blue for segment 8
            else if (i > 97 && i <= 103) colorIndex = 3; // Yellow for segment 9
            else colorIndex = 4; // Dark Gray for segment 10
  
  
            colors.push(...segmentColors[colorIndex]);
        }
    }
  
  
    for (let i = 0; i < segments; i++) {
        for (let j = 0; j < segments; j++) {
            const first = i * (segments + 1) + j;
            const second = first + segments + 1;
  
  
            indices.push(first, second, first + 1);
            indices.push(second, second + 1, first + 1);
        }
    }
  
  
    return { vertices, indices, colors };
  }
  function generateSphere(radius, segments, scale) {
    const vertices = [];
    const indices = [];
    const colors = [];
  
  
    // Define colors for each segment
    const segmentColors = [
        [0, 0, 0],
        [0.8, 0.8, 0.8], // Gray
        [0.5, 0.5, 1.0],
        [1, 1, 0], // Yellow
        [1, 1, 1] // Dark
    ];
  
  
    for (let i = 0; i <= segments; i++) {
        const lat = Math.PI * i / segments;
        const sinLat = Math.sin(lat);
        const cosLat = Math.cos(lat);
  
  
        for (let j = 0; j <= segments; j++) {
            const lng = 2 * Math.PI * j / segments;
            const sinLng = Math.sin(lng);
            const cosLng = Math.cos(lng);
  
  
            const x = cosLng * sinLat * scale;
            const y = cosLat * scale;
            const z = sinLng * sinLat * scale;
  
  
            vertices.push(radius * x, radius * y, radius * z);
  
  
            // Assign color based on segment index
            let colorIndex;
            if (i <= 120) colorIndex = 0; // Blue for segments 1-5
  
  
            else colorIndex = 4; // Dark Gray for segment 10
  
  
            colors.push(...segmentColors[colorIndex]);
        }
    }
  
  
    for (let i = 0; i < segments; i++) {
        for (let j = 0; j < segments; j++) {
            const first = i * (segments + 1) + j;
            const second = first + segments + 1;
  
  
            indices.push(first, second, first + 1);
            indices.push(second, second + 1, first + 1);
        }
    }
  
  
    return { vertices, indices, colors };
  }
  function generateMatahari(radius, segments, scale) {
    const vertices = [];
    const indices = [];
    const colors = [];
  
  
    // Define colors for each segment
    const segmentColors = [
        [1, 206 / 255, 6 / 255]
    ];
  
  
    for (let i = 0; i <= segments; i++) {
        const lat = Math.PI * i / segments;
        const sinLat = Math.sin(lat);
        const cosLat = Math.cos(lat);
  
  
        for (let j = 0; j <= segments; j++) {
            const lng = 2 * Math.PI * j / segments;
            const sinLng = Math.sin(lng);
            const cosLng = Math.cos(lng);
  
  
            const x = cosLng * sinLat * scale;
            const y = cosLat * scale;
            const z = sinLng * sinLat * scale;
  
  
            vertices.push(radius * x, radius * y, radius * z);
  
  
            colors.push(...segmentColors[0]);
        }
    }
  
  
    for (let i = 0; i < segments; i++) {
        for (let j = 0; j < segments; j++) {
            const first = i * (segments + 1) + j;
            const second = first + segments + 1;
  
  
            indices.push(first, second, first + 1);
            indices.push(second, second + 1, first + 1);
        }
    }
  
  
    return { vertices, indices, colors };
  }
  function generateCone(num_triangles, radius, height, color, scale) {
    var cone_vertices = [];
    var cone_colors = [];
    var angle_increment = (2 * Math.PI) / num_triangles;
  
  
    // Generate cone vertices and colors
    for (var i = 0; i < num_triangles; i++) {
        var angle1 = i * angle_increment;
        var angle2 = (i + 1) * angle_increment;
  
  
        // Apex vertex
        cone_vertices.push(0, height, 0);
        cone_colors.push(...color); // White color
  
  
        // First base vertex
        cone_vertices.push(radius * Math.cos(angle1) * scale, 0, radius * Math.sin(angle1) * scale);
        cone_colors.push(...color); // White color
  
  
        // Second base vertex
        cone_vertices.push(radius * Math.cos(angle2) * scale, 0, radius * Math.sin(angle2) * scale);
        cone_colors.push(...color); // White color
    }
  
  
  
  
  
  
    // Define vertices for the bottom of the cone (circle)
    var bottom_center = [0, 0, 0];
    cone_vertices.push(...bottom_center); // Center vertex for bottom
    cone_colors.push(...color); // White color for center vertex
  
  
    for (var i = 0; i < num_triangles; i++) {
        var angle1 = i * angle_increment;
        var angle2 = (i + 1) * angle_increment;
  
  
        cone_vertices.push(radius * Math.cos(angle2) * scale, 0, radius * Math.sin(angle2) * scale); // Vertex on base
        cone_colors.push(...color); // White color for base vertex
  
  
        cone_vertices.push(radius * Math.cos(angle1) * scale, 0, radius * Math.sin(angle1) * scale); // Next vertex on base
        cone_colors.push(...color); // White color for base vertex
  
  
        // Connecting vertex between base and center
        cone_vertices.push(0, 0, 0); // Center vertex
        cone_colors.push(...color); // White color for center vertex
    }
  
  
    return {
        vertices: cone_vertices,
        colors: cone_colors
    };
  }
  function createCylinder(GL, radius, height, num_triangles, r, g, b) {
    var angle_increment = (2 * Math.PI) / num_triangles;
  
  
    // Define cylinder vertices and colors
    var cylinder_vertices = [];
    var top_circle_vertices = [];
    var top_circle_colors = [];
    var bottom_circle_vertices = [];
    var bottom_circle_colors = [];
    var cylinder_colors = [];
  
  
    // Create vertices and colors for cylinder sides
    for (var i = 0; i < num_triangles; i++) {
        var angle1 = i * angle_increment;
        var angle2 = (i + 1) * angle_increment;
  
  
        // First vertex
        cylinder_vertices.push(radius * Math.cos(angle1), 0, radius * Math.sin(angle1));
        cylinder_colors.push(r, g, b);
  
  
        // Second vertex
        cylinder_vertices.push(radius * Math.cos(angle1), height, radius * Math.sin(angle1));
        cylinder_colors.push(r, g, b);
  
  
        // Third vertex
        cylinder_vertices.push(radius * Math.cos(angle2), 0, radius * Math.sin(angle2));
        cylinder_colors.push(r, g, b);
  
  
        // Fourth vertex
        cylinder_vertices.push(radius * Math.cos(angle2), height, radius * Math.sin(angle2));
        cylinder_colors.push(r, g, b);
    }
  
  
    // Create vertices and colors for top and bottom circles
    for (var i = 0; i < num_triangles; i++) {
        var angle1 = i * angle_increment;
        var angle2 = (i + 1) * angle_increment;
  
  
        // Top circle vertices
        top_circle_vertices.push(0, height, 0);
        top_circle_vertices.push(radius * Math.cos(angle1), height, radius * Math.sin(angle1));
        top_circle_vertices.push(radius * Math.cos(angle2), height, radius * Math.sin(angle2));
  
  
        // Top circle colors
        top_circle_colors.push(1, 1, 1);
        top_circle_colors.push(1, 1, 1);
        top_circle_colors.push(1, 1, 1);
  
  
        // Bottom circle vertices
        bottom_circle_vertices.push(0, 0, 0);
        bottom_circle_vertices.push(radius * Math.cos(angle2), 0, radius * Math.sin(angle2));
        bottom_circle_vertices.push(radius * Math.cos(angle1), 0, radius * Math.sin(angle1));
  
  
        // Bottom circle colors
        bottom_circle_colors.push(1, 1, 1);
        bottom_circle_colors.push(1, 1, 1);
        bottom_circle_colors.push(1, 1, 1);
    }
  
  
    // Combine vertices and colors for cylinder, top circle, and bottom circle
    var all_vertices = cylinder_vertices.concat(top_circle_vertices, bottom_circle_vertices);
    var all_colors = cylinder_colors.concat(top_circle_colors, bottom_circle_colors);
  
  
    // Setup VBO for all vertices
    var vbo = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vbo);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(all_vertices), GL.STATIC_DRAW);
  
  
    // Setup VBO for all colors
    var color_vbo = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_vbo);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(all_colors), GL.STATIC_DRAW);
  
  
    return {
        vertices: all_vertices,
        colors: all_colors,
        vbo: vbo,
        color_vbo: color_vbo
    };
  }
  function createCylinderJake(GL, radius, height, num_triangles, r, g, b) {
    var angle_increment = (2 * Math.PI) / num_triangles;
  
  
    // Define cylinder vertices and colors
    var cylinder_vertices = [];
    var top_circle_vertices = [];
    var top_circle_colors = [];
    var bottom_circle_vertices = [];
    var bottom_circle_colors = [];
    var cylinder_colors = [];
  
  
    // Create vertices and colors for cylinder sides
    for (var i = 0; i < num_triangles; i++) {
        var angle1 = i * angle_increment;
        var angle2 = (i + 1) * angle_increment;
  
  
        // First vertex
        cylinder_vertices.push(radius * Math.cos(angle1), 0, radius * Math.sin(angle1));
        cylinder_colors.push(r, g, b);
  
  
        // Second vertex
        cylinder_vertices.push(radius * Math.cos(angle1), height, radius * Math.sin(angle1));
        cylinder_colors.push(r, g, b);
  
  
        // Third vertex
        cylinder_vertices.push(radius * Math.cos(angle2), 0, radius * Math.sin(angle2));
        cylinder_colors.push(r, g, b);
  
  
        // Fourth vertex
        cylinder_vertices.push(radius * Math.cos(angle2), height, radius * Math.sin(angle2));
        cylinder_colors.push(r, g, b);
    }
  
  
    // Create vertices and colors for top and bottom circles
    for (var i = 0; i < num_triangles; i++) {
        var angle1 = i * angle_increment;
        var angle2 = (i + 1) * angle_increment;
  
  
        // Top circle vertices
        top_circle_vertices.push(0, height, 0);
        top_circle_vertices.push(radius * Math.cos(angle1), height, radius * Math.sin(angle1));
        top_circle_vertices.push(radius * Math.cos(angle2), height, radius * Math.sin(angle2));
  
  
        // Top circle colors
        top_circle_colors.push(r, g, b);
        top_circle_colors.push(r, g, b);
        top_circle_colors.push(r, g, b);
  
  
        // Bottom circle vertices
        bottom_circle_vertices.push(0, 0, 0);
        bottom_circle_vertices.push(radius * Math.cos(angle2), 0, radius * Math.sin(angle2));
        bottom_circle_vertices.push(radius * Math.cos(angle1), 0, radius * Math.sin(angle1));
  
  
        // Bottom circle colors
        bottom_circle_colors.push(r, g, b);
        bottom_circle_colors.push(r, g, b);
        bottom_circle_colors.push(r, g, b);
    }
  
  
    // Combine vertices and colors for cylinder, top circle, and bottom circle
    var all_vertices = cylinder_vertices.concat(top_circle_vertices, bottom_circle_vertices);
    var all_colors = cylinder_colors.concat(top_circle_colors, bottom_circle_colors);
  
  
    // Setup VBO for all vertices
    var vbo = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, vbo);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(all_vertices), GL.STATIC_DRAW);
  
  
    // Setup VBO for all colors
    var color_vbo = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, color_vbo);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(all_colors), GL.STATIC_DRAW);
  
  
    return {
        vertices: all_vertices,
        colors: all_colors,
        vbo: vbo,
        color_vbo: color_vbo
    };
  }
  function createKuboid(GL, width, height, depth) {
    var vertices = [
        // Front face
        -width, -height, depth,
        width, -height, depth,
        width, height, depth,
        -width, height, depth,
  
  
        // Back face
        -width, -height, -depth,
        -width, height, -depth,
        width, height, -depth,
        width, -height, -depth,
  
  
        // Top face
        -width, height, -depth,
        -width, height, depth,
        width, height, depth,
        width, height, -depth,
  
  
        // Bottom face
        -width, -height, -depth,
        width, -height, -depth,
        width, -height, depth,
        -width, -height, depth,
  
  
        // Right face
        width, -height, -depth,
        width, height, -depth,
        width, height, depth,
        width, -height, depth,
  
  
        // Left face
        -width, -height, -depth,
        -width, -height, depth,
        -width, height, depth,
        -width, height, -depth
    ];
  
  
    var colors = [
        // Front face
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
  
  
        // Back face
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
  
  
        // Top face
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
  
  
        // Bottom face
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
  
  
        // Right face
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
  
  
        // Left face
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
    ];
  
  
    var indices = [
        0, 1, 2, 0, 2, 3, // Front face
        4, 5, 6, 4, 6, 7, // Back face
        8, 9, 10, 8, 10, 11, // Top face
        12, 13, 14, 12, 14, 15, // Bottom face
        16, 17, 18, 16, 18, 19, // Right face
        20, 21, 22, 20, 22, 23 // Left face
    ];
  
  
    return {
        vertices: vertices,
        colors: colors,
        indices: indices
    };
  }
  function createLemari(GL, width, height, depth) {
    var vertices = [
        // Back face
        -width, -height, -depth,
        -width, height, -depth,
        width, height, -depth,
        width, -height, -depth,
  
  
        // Top face
        -width, height, -depth,
        -width, height, depth,
        width, height, depth,
        width, height, -depth,
  
  
        // Bottom face
        -width, -height, -depth,
        width, -height, -depth,
        width, -height, depth,
        -width, -height, depth,
  
  
        // Right face
        width, -height, -depth,
        width, height, -depth,
        width, height, depth,
        width, -height, depth,
  
  
        // Left face
        -width, -height, -depth,
        -width, -height, depth,
        -width, height, depth,
        -width, height, -depth
    ];
  
  
    var colors = [
        // Back face
        0.6, 0.3, 0.0,
        0.6, 0.3, 0.0,
        0.6, 0.3, 0.0,
        0.6, 0.3, 0.0,
  
  
        // Top face
        0.6, 0.35, 0.0,
        0.6, 0.35, 0.0,
        0.6, 0.35, 0.0,
        0.6, 0.35, 0.0,
  
  
        // Bottom face
        0.6, 0.3, 0.0,
        0.6, 0.3, 0.0,
        0.6, 0.3, 0.0,
        0.6, 0.3, 0.0,
  
  
        // Right face
        0.65, 0.35, 0.0,
        0.65, 0.35, 0.0,
        0.65, 0.35, 0.0,
        0.65, 0.35, 0.0,
  
  
        // Left face
        0.55, 0.35, 0.0,
        0.55, 0.35, 0.0,
        0.55, 0.35, 0.0,
        0.55, 0.35, 0.0,
    ];
  
  
    var indices = [
        //depannya bolong
        1, 2, 3, 1, 3, 0, // Back face
        5, 6, 7, 5, 7, 4, // Top face
        9, 10, 11, 9, 11, 8, // Bottom face
        13, 14, 15, 13, 15, 12, // Right face
        17, 18, 19, 17, 19, 16 // Left face
  
  
    ];
  
  
    return {
        vertices: vertices,
        colors: colors,
        indices: indices
    };
  }
  
  
  
  
  function main() {
    var CANVAS = document.getElementById("myCanvas");
    CANVAS.width = window.innerWidth;
    CANVAS.height = window.innerHeight;
    var GL;
    try {
        GL = CANVAS.getContext("webgl", { antialias: true });
    } catch (e) {
        alert("WebGL context cannot be initialized");
        return false;
    }
    //shaders
    var shader_vertex_source = `
      attribute vec3 position;
      attribute vec3 color; // New attribute for color
   
      uniform mat4 PMatrix;
      uniform mat4 VMatrix;
      uniform mat4 MMatrix;
   
      varying vec3 vColor; // Varying variable to pass color to fragment shader
   
      void main(void) {
          gl_Position = PMatrix * VMatrix * MMatrix * vec4(position, 1.0);
          vColor = color; // Pass color to fragment shader
      }
    `;
    var shader_fragment_source = `
        precision mediump float;
        varying vec3 vColor; // Receive color from the vertex shader
        void main(void) {
            gl_FragColor = vec4(vColor, 1.0); // Use received color for fragment
        }
    `;
    var compile_shader = function (source, type, typeString) {
        var shader = GL.createShader(type);
        GL.shaderSource(shader, source);
        GL.compileShader(shader);
        if (!GL.getShaderParameter(shader, GL.COMPILE_STATUS)) {
            alert("ERROR IN " + typeString + " SHADER: " + GL.getShaderInfoLog(shader));
            return false;
        }
        return shader;
    };
  
  
    var shader_vertex = compile_shader(shader_vertex_source, GL.VERTEX_SHADER, "VERTEX");
    var shader_fragment = compile_shader(shader_fragment_source, GL.FRAGMENT_SHADER, "FRAGMENT");
    var SHADER_PROGRAM = GL.createProgram();
    GL.attachShader(SHADER_PROGRAM, shader_vertex);
    GL.attachShader(SHADER_PROGRAM, shader_fragment);
    GL.linkProgram(SHADER_PROGRAM);
  
  
  
  
    var _color = GL.getAttribLocation(SHADER_PROGRAM, "color");
    var _position = GL.getAttribLocation(SHADER_PROGRAM, "position");
    //uniform
    var _PMatrix = GL.getUniformLocation(SHADER_PROGRAM, "PMatrix"); //projection
    var _VMatrix = GL.getUniformLocation(SHADER_PROGRAM, "VMatrix"); //View
    var _MMatrix = GL.getUniformLocation(SHADER_PROGRAM, "MMatrix"); //Model
  
  
    GL.enableVertexAttribArray(_color);
    GL.enableVertexAttribArray(_position);
    GL.useProgram(SHADER_PROGRAM);
  
  
  
  
  
  
  
  
    //##############################################################################################################
    //########################################### OBJECT'S #########################################################
    //##############################################################################################################
  
  
    //===================================== POHON 1 ============================================================
    //------------------------------------ BATANG 1 ------------------------------------------------------------
    var batangpohon1 = createCylinder(GL, 0.7, 10, 50, 150 / 255, 75 / 255, 0);
    console.log(batangpohon1.vertices);
  
  
    var BATANG_POHON1_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BATANG_POHON1_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(batangpohon1.vertices), GL.STATIC_DRAW);
  
  
    var BATANG_POHON1_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BATANG_POHON1_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(batangpohon1.colors), GL.STATIC_DRAW);
  
  
    var MODEL_BATANG_POHON1_MATRIX = LIBS.get_I4();
    LIBS.translateX(MODEL_BATANG_POHON1_MATRIX, LIBS.degToRad(-860));
    LIBS.translateY(MODEL_BATANG_POHON1_MATRIX, LIBS.degToRad(120));
    LIBS.translateZ(MODEL_BATANG_POHON1_MATRIX, LIBS.degToRad(-500)); // ATUR TINGGI
    LIBS.rotateX(MODEL_BATANG_POHON1_MATRIX, LIBS.degToRad(90));
    //------------------------------------ BATANG 1 ------------------------------------------------------------
    //------------------------------------- DAUN 1 -------------------------------------------------------------
    var daunpohon1 = generateElipticPara(5, 100, 1, 11 / 255, 102 / 255, 35 / 255, true);
  
  
    var DAUN_POHON1_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, DAUN_POHON1_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(daunpohon1.vertices), GL.STATIC_DRAW);
  
  
    var DAUN_POHON1_FACES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, DAUN_POHON1_FACES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(daunpohon1.indices), GL.STATIC_DRAW);
  
  
    var DAUN_POHON1_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, DAUN_POHON1_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(daunpohon1.colors), GL.STATIC_DRAW);
  
  
    var MODEL_DAUN_POHON1_MATRIX = LIBS.get_I4();
    LIBS.translateX(MODEL_DAUN_POHON1_MATRIX, LIBS.degToRad(-860));
    LIBS.translateY(MODEL_DAUN_POHON1_MATRIX, LIBS.degToRad(120));
    LIBS.translateZ(MODEL_DAUN_POHON1_MATRIX, LIBS.degToRad(-1000));
    LIBS.rotateX(MODEL_DAUN_POHON1_MATRIX, LIBS.degToRad(0));
    //---------------------------------------------------------------------------------------------------------
    //==========================================================================================================//
  
  
  
  
    //============================================== KALENG ====================================================//
    //--------------------------------------------- KALENG 1 ---------------------------------------------------- 
    var kaleng1 = createCylinder(GL, 0.2, 1, 50, 1, 0, 0);
    console.log(kaleng1.vertices);
  
  
    var KALENG1_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, KALENG1_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(kaleng1.vertices), GL.STATIC_DRAW);
  
  
    var KALENG1_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, KALENG1_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(kaleng1.colors), GL.STATIC_DRAW);
  
  
    // LOKASI KALENG1: 
    var MODEL_KALENG1_MATRIX = LIBS.get_I4();
    LIBS.translateX(MODEL_KALENG1_MATRIX, LIBS.degToRad(-230));
    LIBS.translateY(MODEL_KALENG1_MATRIX, LIBS.degToRad(-500));
    LIBS.translateZ(MODEL_KALENG1_MATRIX, LIBS.degToRad(-583));
    LIBS.rotateX(MODEL_KALENG1_MATRIX, LIBS.degToRad(90));
    //----------------------------------------------------------------------------------------------------------- 
  
  
  
  
    //--------------------------------------------- KALENG 2 ---------------------------------------------------- 
    var kaleng2 = createCylinder(GL, 0.2, 1, 50, 1, 0, 0);
    console.log(kaleng2.vertices);
  
  
    var KALENG2_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, KALENG2_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(kaleng2.vertices), GL.STATIC_DRAW);
  
  
    var KALENG2_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, KALENG2_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(kaleng2.colors), GL.STATIC_DRAW);
  
  
    // LOKASI KALENG2: 
    var MODEL_KALENG2_MATRIX = LIBS.get_I4();
    LIBS.translateX(MODEL_KALENG2_MATRIX, LIBS.degToRad(-230));
    LIBS.translateY(MODEL_KALENG2_MATRIX, LIBS.degToRad(-450));
    LIBS.translateZ(MODEL_KALENG2_MATRIX, LIBS.degToRad(-583));
    LIBS.rotateX(MODEL_KALENG2_MATRIX, LIBS.degToRad(90));
    //----------------------------------------------------------------------------------------------------------- 
    //==========================================================================================================//
  
  
  
  
  
  
  
  
    //==================================== LEMARI ========================================
    // Setup VAO
    var position_vao = GL.getAttribLocation(SHADER_PROGRAM, "position");
    var color_vao = GL.getAttribLocation(SHADER_PROGRAM, "color");
    //------------------------------ Create kiri tangga ----------------------------------
    var kiri_tangga = createKuboid(GL, 0.05, 6, 0.05);
    var kiri_tangga_vertex = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, kiri_tangga_vertex);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(kiri_tangga.vertices), GL.STATIC_DRAW);
  
  
    var kiri_tangga_color = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, kiri_tangga_color);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(kiri_tangga.colors), GL.STATIC_DRAW);
  
  
    var kiri_tangga_faces = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, kiri_tangga_faces);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(kiri_tangga.indices), GL.STATIC_DRAW);
  
  
    var kanan_tangga = createKuboid(GL, 0.05, 6, 0.05);
    var kanan_tangga_vertex = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, kanan_tangga_vertex);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(kanan_tangga.vertices), GL.STATIC_DRAW);
  
  
    var kanan_tangga_color = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, kanan_tangga_color);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(kanan_tangga.colors), GL.STATIC_DRAW);
  
  
    var kanan_tangga_faces = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, kanan_tangga_faces);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(kanan_tangga.indices), GL.STATIC_DRAW);
    //----------------------------------------------------------------------------------------
  
  
  
  
    //----------------------------- Create tengah tangga -------------------------------------
    var tengah_tangga_1 = createKuboid(GL, 0.05, 0.8, 0.05);
    var tengah_tangga_1_vertex = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_1_vertex);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tengah_tangga_1.vertices), GL.STATIC_DRAW);
  
  
    var tengah_tangga_1_color = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_1_color);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tengah_tangga_1.colors), GL.STATIC_DRAW);
  
  
    var tengah_tangga_1_faces = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, tengah_tangga_1_faces);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(tengah_tangga_1.indices), GL.STATIC_DRAW);
  
  
    var tengah_tangga_2 = createKuboid(GL, 0.05, 0.8, 0.05);
    var tengah_tangga_2_vertex = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_2_vertex);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tengah_tangga_2.vertices), GL.STATIC_DRAW);
  
  
    var tengah_tangga_2_color = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_2_color);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tengah_tangga_2.colors), GL.STATIC_DRAW);
  
  
    var tengah_tangga_2_faces = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, tengah_tangga_2_faces);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(tengah_tangga_2.indices), GL.STATIC_DRAW);
  
  
    var tengah_tangga_3 = createKuboid(GL, 0.05, 0.8, 0.05);
    var tengah_tangga_4 = createKuboid(GL, 0.05, 0.8, 0.05);
    var tengah_tangga_5 = createKuboid(GL, 0.05, 0.8, 0.05);
    var tengah_tangga_6 = createKuboid(GL, 0.05, 0.8, 0.05);
    var tengah_tangga_7 = createKuboid(GL, 0.05, 0.8, 0.05);
    var tengah_tangga_8 = createKuboid(GL, 0.05, 0.8, 0.05);
    var tengah_tangga_9 = createKuboid(GL, 0.05, 0.8, 0.05);
    var tengah_tangga_10 = createKuboid(GL, 0.05, 0.8, 0.05);
  
  
    var tengah_tangga_3_vertex = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_3_vertex);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tengah_tangga_3.vertices), GL.STATIC_DRAW);
  
  
    var tengah_tangga_3_color = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_3_color);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tengah_tangga_3.colors), GL.STATIC_DRAW);
  
  
    var tengah_tangga_3_faces = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, tengah_tangga_3_faces);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(tengah_tangga_3.indices), GL.STATIC_DRAW);
  
  
    var tengah_tangga_4_vertex = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_4_vertex);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tengah_tangga_4.vertices), GL.STATIC_DRAW);
  
  
    var tengah_tangga_4_color = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_4_color);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tengah_tangga_4.colors), GL.STATIC_DRAW);
  
  
    var tengah_tangga_4_faces = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, tengah_tangga_4_faces);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(tengah_tangga_4.indices), GL.STATIC_DRAW);
  
  
    var tengah_tangga_5_vertex = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_5_vertex);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tengah_tangga_5.vertices), GL.STATIC_DRAW);
  
  
    var tengah_tangga_5_color = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_5_color);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tengah_tangga_5.colors), GL.STATIC_DRAW);
  
  
    var tengah_tangga_5_faces = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, tengah_tangga_5_faces);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(tengah_tangga_5.indices), GL.STATIC_DRAW);
  
  
    var tengah_tangga_6_vertex = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_6_vertex);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tengah_tangga_6.vertices), GL.STATIC_DRAW);
  
  
    var tengah_tangga_6_color = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_6_color);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tengah_tangga_6.colors), GL.STATIC_DRAW);
  
  
    var tengah_tangga_6_faces = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, tengah_tangga_6_faces);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(tengah_tangga_6.indices), GL.STATIC_DRAW);
  
  
    var tengah_tangga_7_vertex = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_7_vertex);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tengah_tangga_7.vertices), GL.STATIC_DRAW);
  
  
    var tengah_tangga_7_color = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_7_color);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tengah_tangga_7.colors), GL.STATIC_DRAW);
  
  
    var tengah_tangga_7_faces = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, tengah_tangga_7_faces);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(tengah_tangga_7.indices), GL.STATIC_DRAW);
  
  
    var tengah_tangga_8_vertex = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_8_vertex);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tengah_tangga_8.vertices), GL.STATIC_DRAW);
  
  
    var tengah_tangga_8_color = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_8_color);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tengah_tangga_8.colors), GL.STATIC_DRAW);
  
  
    var tengah_tangga_8_faces = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, tengah_tangga_8_faces);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(tengah_tangga_8.indices), GL.STATIC_DRAW);
  
  
    var tengah_tangga_9_vertex = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_9_vertex);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tengah_tangga_9.vertices), GL.STATIC_DRAW);
  
  
    var tengah_tangga_9_color = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_9_color);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tengah_tangga_9.colors), GL.STATIC_DRAW);
  
  
    var tengah_tangga_9_faces = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, tengah_tangga_9_faces);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(tengah_tangga_9.indices), GL.STATIC_DRAW);
  
  
    var tengah_tangga_10_vertex = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_10_vertex);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tengah_tangga_10.vertices), GL.STATIC_DRAW);
  
  
    var tengah_tangga_10_color = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_10_color);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tengah_tangga_10.colors), GL.STATIC_DRAW);
  
  
    var tengah_tangga_10_faces = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, tengah_tangga_10_faces);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(tengah_tangga_10.indices), GL.STATIC_DRAW);
    //----------------------------------------------------------------------------------------
  
  
    //----------------------------------- Create lemari --------------------------------------
    var lemari = createLemari(GL, 3, 5, 0.8);
    var lemari_vertex = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, lemari_vertex);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(lemari.vertices), GL.STATIC_DRAW);
  
  
    var lemari_color = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, lemari_color);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(lemari.colors), GL.STATIC_DRAW);
  
  
    var lemari_faces = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, lemari_faces);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(lemari.indices), GL.STATIC_DRAW);
  
  
    var laci_lemari_1 = createKuboid(GL, 2.86, 0.05, 0.8);
    var laci_lemari_1_vertex = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, laci_lemari_1_vertex);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(laci_lemari_1.vertices), GL.STATIC_DRAW);
  
  
    var laci_lemari_1_color = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, laci_lemari_1_color);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(laci_lemari_1.colors), GL.STATIC_DRAW);
  
  
    var laci_lemari_1_faces = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, laci_lemari_1_faces);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(laci_lemari_1.indices), GL.STATIC_DRAW);
  
  
    var laci_lemari_2 = createKuboid(GL, 2.9, 0.05, 0.8);
    var laci_lemari_2_vertex = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, laci_lemari_2_vertex);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(laci_lemari_2.vertices), GL.STATIC_DRAW);
  
  
    var laci_lemari_2_color = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, laci_lemari_2_color);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(laci_lemari_2.colors), GL.STATIC_DRAW);
  
  
    var laci_lemari_2_faces = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, laci_lemari_2_faces);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(laci_lemari_2.indices), GL.STATIC_DRAW);
    //----------------------------------------------------------------------------------------
  
  
  
  
    //---------------------------------- Matrix ----------------------------------------------
    var PROJECTION_MATRIX = LIBS.get_projection(40, CANVAS.width / CANVAS.height, 1, 100);
    var VIEW_MATRIX = LIBS.get_I4();
    var MODEL_MATRIX = LIBS.get_I4();
  
  
  
  
    var kiri_tangga_MODEL_MATRIX = LIBS.get_I4();
  
  
    var tengah_tangga_1_MODEL_MATRIX = LIBS.get_I4();
    var tengah_tangga_2_MODEL_MATRIX = LIBS.get_I4();
    var tengah_tangga_3_MODEL_MATRIX = LIBS.get_I4();
    var tengah_tangga_4_MODEL_MATRIX = LIBS.get_I4();
    var tengah_tangga_5_MODEL_MATRIX = LIBS.get_I4();
    var tengah_tangga_6_MODEL_MATRIX = LIBS.get_I4();
    var tengah_tangga_7_MODEL_MATRIX = LIBS.get_I4();
    var tengah_tangga_8_MODEL_MATRIX = LIBS.get_I4();
    var tengah_tangga_9_MODEL_MATRIX = LIBS.get_I4();
    var tengah_tangga_10_MODEL_MATRIX = LIBS.get_I4();
  
  
    var kanan_tangga_MODEL_MATRIX = LIBS.get_I4();
  
  
    var lemari_MODEL_MATRIX = LIBS.get_I4();
    var laci_lemari_1_MODEL_MATRIX = LIBS.get_I4();
    var laci_lemari_2_MODEL_MATRIX = LIBS.get_I4();
    //----------------------------------------------------------------------------------------
  
  
  
  
    //-------------------------- Atur tengah tangga ------------------------------------------
    LIBS.rotateZ(tengah_tangga_1_MODEL_MATRIX, LIBS.degToRad(90));
    LIBS.translateX(tengah_tangga_1_MODEL_MATRIX, 0.8);
    LIBS.translateY(tengah_tangga_1_MODEL_MATRIX, 5);
  
  
    LIBS.rotateZ(tengah_tangga_2_MODEL_MATRIX, LIBS.degToRad(90));
    LIBS.translateX(tengah_tangga_2_MODEL_MATRIX, 0.8);
    LIBS.translateY(tengah_tangga_2_MODEL_MATRIX, 4);
  
  
    LIBS.rotateZ(tengah_tangga_3_MODEL_MATRIX, LIBS.degToRad(90));
    LIBS.translateX(tengah_tangga_3_MODEL_MATRIX, 0.8);
    LIBS.translateY(tengah_tangga_3_MODEL_MATRIX, 3);
  
  
    LIBS.rotateZ(tengah_tangga_4_MODEL_MATRIX, LIBS.degToRad(90));
    LIBS.translateX(tengah_tangga_4_MODEL_MATRIX, 0.8);
    LIBS.translateY(tengah_tangga_4_MODEL_MATRIX, 2);
  
  
    LIBS.rotateZ(tengah_tangga_5_MODEL_MATRIX, LIBS.degToRad(90));
    LIBS.translateX(tengah_tangga_5_MODEL_MATRIX, 0.8);
    LIBS.translateY(tengah_tangga_5_MODEL_MATRIX, 1);
  
  
    LIBS.rotateZ(tengah_tangga_6_MODEL_MATRIX, LIBS.degToRad(90));
    LIBS.translateX(tengah_tangga_6_MODEL_MATRIX, 0.8);
    LIBS.translateY(tengah_tangga_6_MODEL_MATRIX, 0);
  
  
    LIBS.rotateZ(tengah_tangga_7_MODEL_MATRIX, LIBS.degToRad(90));
    LIBS.translateX(tengah_tangga_7_MODEL_MATRIX, 0.8);
    LIBS.translateY(tengah_tangga_7_MODEL_MATRIX, -1.4);
  
  
    LIBS.rotateZ(tengah_tangga_8_MODEL_MATRIX, LIBS.degToRad(90));
    LIBS.translateX(tengah_tangga_8_MODEL_MATRIX, 0.8);
    LIBS.translateY(tengah_tangga_8_MODEL_MATRIX, -2.6);
  
  
    LIBS.rotateZ(tengah_tangga_9_MODEL_MATRIX, LIBS.degToRad(90));
    LIBS.translateX(tengah_tangga_9_MODEL_MATRIX, 0.8);
    LIBS.translateY(tengah_tangga_9_MODEL_MATRIX, -3.8);
  
  
    LIBS.rotateZ(tengah_tangga_10_MODEL_MATRIX, LIBS.degToRad(90));
    LIBS.translateX(tengah_tangga_10_MODEL_MATRIX, 0.8);
    LIBS.translateY(tengah_tangga_10_MODEL_MATRIX, -5);
    //----------------------------------------------------------------------------------------
  
  
    //---------------------- atur kanan tangga -----------------------------------------------
    LIBS.translateX(kanan_tangga_MODEL_MATRIX, 1.6);
    //----------------------------------------------------------------------------------------
  
  
    //---------------------- atur laci lemari -----------------------------------------------
    LIBS.translateY(laci_lemari_1_MODEL_MATRIX, 2);
    LIBS.translateY(laci_lemari_2_MODEL_MATRIX, -2);
    //----------------------------------------------------------------------------------------
    //========================================================================================
  
  
  
  
  
  
    //==================================== MATAHARI =======================================================//
    var matahari = generateMatahari(1, 300, 10);
    console.log(matahari.vertices);
    console.log(matahari.indices);
  
  
    var MATAHARI_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, MATAHARI_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(matahari.vertices), GL.STATIC_DRAW);
  
  
    var MATAHARI_FACES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, MATAHARI_FACES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(matahari.indices), GL.STATIC_DRAW);
  
  
    var MATAHARI_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, MATAHARI_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(matahari.colors), GL.STATIC_DRAW);
  
  
    var MODEL_MATAHARI_MATRIX = LIBS.get_I4();
    LIBS.translateX(MODEL_MATAHARI_MATRIX, LIBS.degToRad(900));
    LIBS.translateY(MODEL_MATAHARI_MATRIX, LIBS.degToRad(900));
    LIBS.translateZ(MODEL_MATAHARI_MATRIX, LIBS.degToRad(-2000));
    //=========================================================================================================//
  
  
  
  
  
  
    //========================== KINCIR ANGIN ================================================================//
    //----------------------- Bangunan Kincir Bawah ----------------------------------------------------------// 
    var kincirbawah = generateElipticPara(5, 200, 1, 254/255, 230/255, 168/255);
  
  
    var KINCIR_BAWAH_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, KINCIR_BAWAH_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(kincirbawah.vertices), GL.STATIC_DRAW);
  
  
    var KINCIR_BAWAH_FACES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, KINCIR_BAWAH_FACES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(kincirbawah.indices), GL.STATIC_DRAW);
  
  
    var KINCIR_BAWAH_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, KINCIR_BAWAH_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(kincirbawah.colors), GL.STATIC_DRAW);
  
  
    var MODEL_KINCIR_BAWAH_MATRIX = LIBS.get_I4();
    LIBS.translateX(MODEL_KINCIR_BAWAH_MATRIX, LIBS.degToRad(-2000));
    LIBS.translateY(MODEL_KINCIR_BAWAH_MATRIX, LIBS.degToRad(-1100));
    LIBS.translateZ(MODEL_KINCIR_BAWAH_MATRIX, LIBS.degToRad(-600));
    LIBS.rotateX(MODEL_KINCIR_BAWAH_MATRIX, LIBS.degToRad(0));
    //--------------------------------------------------------------------------------------------------------// 
  
  
  
  
    //----------------------- Bangunan Kincir Atas ----------------------------------------------------------// 
    var kinciratas = generateElipticPara(5, 270, 1.5, 178/255, 34/255, 3/255);
  
  
    var KINCIR_ATAS_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, KINCIR_ATAS_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(kinciratas.vertices), GL.STATIC_DRAW);
  
  
    var KINCIR_ATAS_FACES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, KINCIR_ATAS_FACES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(kinciratas.indices), GL.STATIC_DRAW);
  
  
    var KINCIR_ATAS_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, KINCIR_ATAS_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(kinciratas.colors), GL.STATIC_DRAW);
  
  
    var MODEL_KINCIR_ATAS_MATRIX = LIBS.get_I4();
    LIBS.translateX(MODEL_KINCIR_ATAS_MATRIX, LIBS.degToRad(-2000));
    LIBS.translateY(MODEL_KINCIR_ATAS_MATRIX, LIBS.degToRad(-1100));
    LIBS.translateZ(MODEL_KINCIR_ATAS_MATRIX, LIBS.degToRad(-1200));
    LIBS.rotateX(MODEL_KINCIR_ATAS_MATRIX, LIBS.degToRad(0));
    //--------------------------------------------------------------------------------------------------------// 
  
  
  
  
    //--------------------------------------------- TABUNG KINCIR ---------------------------------------------------- 
    var tabungkincir = createCylinder(GL, 1, 5, 50, 1, 0, 0);
    console.log(tabungkincir.vertices);
  
  
    var TABUNG_KINCIR_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, TABUNG_KINCIR_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tabungkincir.vertices), GL.STATIC_DRAW);
  
  
    var TABUNG_KINCIR_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, TABUNG_KINCIR_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tabungkincir.colors), GL.STATIC_DRAW);
  
  
    // LOKASI KINCIR: 
    var MODEL_TABUNG_KINCIR_MATRIX = LIBS.get_I4();
    LIBS.translateX(MODEL_TABUNG_KINCIR_MATRIX, LIBS.degToRad(-1920));
    LIBS.translateY(MODEL_TABUNG_KINCIR_MATRIX, LIBS.degToRad(-900));
    LIBS.translateZ(MODEL_TABUNG_KINCIR_MATRIX, LIBS.degToRad(-800));
    LIBS.rotateX(MODEL_TABUNG_KINCIR_MATRIX, LIBS.degToRad(0));
    LIBS.rotateZ(MODEL_TABUNG_KINCIR_MATRIX, LIBS.degToRad(-40));
    //----------------------------------------------------------------------------------------------------------- 
  
  
    //--------------------------------------------- KINCIR 1 ------------------------------------------------//
    var kincir1 = createCylinder(GL, 0.5, 7, 50, 0, 0, 0);
    console.log(kincir1.vertices);
  
  
    var KINCIR1_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, KINCIR1_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(kincir1.vertices), GL.STATIC_DRAW);
  
  
    var KINCIR1_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, KINCIR1_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(kincir1.colors), GL.STATIC_DRAW);
  
  
    var MODEL_KINCIR1_MATRIX = LIBS.get_I4();
    LIBS.translateX(MODEL_KINCIR1_MATRIX, LIBS.degToRad(-1780)); // Move the cylinder along the X-axis
    LIBS.translateY(MODEL_KINCIR1_MATRIX, LIBS.degToRad(-750)); // Move the cylinder along the Y-axis
    LIBS.translateZ(MODEL_KINCIR1_MATRIX, LIBS.degToRad(-800)); // Move the cylinder along the Z-axis
    LIBS.rotateZ(MODEL_KINCIR1_MATRIX, -1.5);
    LIBS.rotateZ(MODEL_KINCIR1_MATRIX, LIBS.degToRad(-40));
    //------------------------------------------------------------------------------------------------------/
    
    //--------------------------------------------- KINCIR 2 ------------------------------------------------//
    var kincir2 = createCylinder(GL, 0.5, 7, 50, 0, 0, 0);
    console.log(kincir2.vertices);
  
  
    var KINCIR2_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, KINCIR2_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(kincir2.vertices), GL.STATIC_DRAW);
  
  
    var KINCIR2_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, KINCIR2_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(kincir2.colors), GL.STATIC_DRAW);
  
  
    var MODEL_KINCIR2_MATRIX = LIBS.get_I4();
    LIBS.translateX(MODEL_KINCIR2_MATRIX, LIBS.degToRad(-1780)); // Move the cylinder along the X-axis
    LIBS.translateY(MODEL_KINCIR2_MATRIX, LIBS.degToRad(-750)); // Move the cylinder along the Y-axis
    LIBS.translateZ(MODEL_KINCIR2_MATRIX, LIBS.degToRad(-800)); // Move the cylinder along the Z-axis
    LIBS.rotateZ(MODEL_KINCIR2_MATRIX, -1.5);
    LIBS.rotateY(MODEL_KINCIR2_MATRIX, 2);
    LIBS.rotateZ(MODEL_KINCIR2_MATRIX, LIBS.degToRad(-40));
  
  
    //------------------------------------------------------------------------------------------------------/
  
  
    //--------------------------------------------- KINCIR 3 ------------------------------------------------//
    var kincir3 = createCylinder(GL, 0.5, 7, 50, 0, 0, 0);
    console.log(kincir3.vertices);
  
  
    var KINCIR3_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, KINCIR3_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(kincir3.vertices), GL.STATIC_DRAW);
  
  
    var KINCIR3_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, KINCIR3_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(kincir3.colors), GL.STATIC_DRAW);
  
  
    var MODEL_KINCIR3_MATRIX = LIBS.get_I4();
    LIBS.translateX(MODEL_KINCIR3_MATRIX, LIBS.degToRad(-1780)); // Move the cylinder along the X-axis
    LIBS.translateY(MODEL_KINCIR3_MATRIX, LIBS.degToRad(-750)); // Move the cylinder along the Y-axis
    LIBS.translateZ(MODEL_KINCIR3_MATRIX, LIBS.degToRad(-800)); // Move the cylinder along the Z-axis
    LIBS.rotateZ(MODEL_KINCIR3_MATRIX, -1.5);
    LIBS.rotateY(MODEL_KINCIR3_MATRIX, 4);
    LIBS.rotateZ(MODEL_KINCIR3_MATRIX, LIBS.degToRad(-40));
  
  
    //------------------------------------------------------------------------------------------------------/
  
  
  
  
  
  
  
  
    //##############################################################################################################
    //########################################### OBJECT'S #########################################################
    //##############################################################################################################
  
  
  
  
  
  
    //##############################################################################################################
    //################################# TV The Son of Jake #########################################################
    //##############################################################################################################
    /*========================= THE SPHERE (TV's BODY) ========================= */
    var sphere = createSphere(1, 300);
    console.log(sphere.vertices);
    console.log(sphere.indices);
  
  
    var TRIANGLE_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, TRIANGLE_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(sphere.vertices), GL.STATIC_DRAW);
  
  
    var TRIANGLE_FACES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, TRIANGLE_FACES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(sphere.indices), GL.STATIC_DRAW);
  
  
    var TRIANGLE_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, TRIANGLE_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(sphere.colors), GL.STATIC_DRAW);
    /*================================================================ */
  
  
  
  
  
  
    /*========================= THE CONE ========================= */
    var cone = generateCone(300, 0.5, 1.4, [1.0, 1.0, 1.0], 1); // Adjust parameters as needed
  
  
    var CONE_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, CONE_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(cone.vertices), GL.STATIC_DRAW);
  
  
    var CONE_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, CONE_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(cone.colors), GL.STATIC_DRAW);
    /*================================================================ */
  
  
  
  
  
  
    /*========================= THE Left Eye ========================= */
    var eye = createCylinder(GL, 0.11, -0.2, 50, 0, 0, 0);
    console.log(eye.vertices);
  
  
  
  
    var EYE_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, EYE_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(eye.vertices), GL.STATIC_DRAW);
  
  
    var EYE_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, EYE_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(eye.colors), GL.STATIC_DRAW);
    /*================================================================ */
  
  
  
  
  
  
    /*========================= The Right Eye ===============================*/
    var second_eye = createCylinder(GL, 0.11, -0.2, 50, 0, 0, 0);
    console.log(second_eye.vertices);
  
  
    var SECOND_EYE_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, SECOND_EYE_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(second_eye.vertices), GL.STATIC_DRAW);
  
  
    var SECOND_EYE_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, SECOND_EYE_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(second_eye.colors), GL.STATIC_DRAW);
    /*================================================================ */
  
  
  
  
  
  
    /*========================= Alis kanan ===============================*/
    var right_brow = createKuboid(GL, 0.13, 0.03, 0.02);
    console.log(right_brow.vertices);
  
  
    var RIGHT_BROW_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_BROW_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(right_brow.vertices), GL.STATIC_DRAW);
  
  
    var RIGHT_BROW_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_BROW_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(right_brow.colors), GL.STATIC_DRAW);
    /*================================================================ */
  
  
  
  
  
  
    /*========================= Alis kiri ===============================*/
    var left_brow = createKuboid(GL, 0.13, 0.03, 0.02);
    console.log(left_brow.vertices);
  
  
    var LEFT_BROW_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_BROW_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(left_brow.vertices), GL.STATIC_DRAW);
  
  
    var LEFT_BROW_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_BROW_COLOR);
  
  
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(left_brow.colors), GL.STATIC_DRAW);
    /*================================================================ */
  
  
  
  
  
  
    /*========================= Tangan kiri ===============================*/
    var left_hand = generateEllipsoid(0, 0, 0, 0.06, 0.3, 0.06, 50, [0.5, 0.5, 1.0]);
    console.log(left_hand.vertices); // Check the generated vertices
  
  
    var LEFT_HAND_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_HAND_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(left_hand.vertices), GL.STATIC_DRAW);
  
  
    var LEFT_HAND_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, LEFT_HAND_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(left_hand.indices), GL.STATIC_DRAW);
  
  
    var LEFT_HAND_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_HAND_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(left_hand.colors), GL.STATIC_DRAW);
    /*====================================================================== */
  
  
  
  
  
  
    /*========================= Tangan kanan ===============================*/
    var right_hand = generateEllipsoid(0, 0, 0, 0.06, 0.3, 0.06, 50, [0.5, 0.5, 1.0]);
    console.log(right_hand.vertices); // Check the generated vertices
  
  
    var RIGHT_HAND_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_HAND_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(right_hand.vertices), GL.STATIC_DRAW);
  
  
    var RIGHT_HAND_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, RIGHT_HAND_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(right_hand.indices), GL.STATIC_DRAW);
  
  
    var RIGHT_HAND_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_HAND_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(right_hand.colors), GL.STATIC_DRAW);
    /*================================================================ */
  
  
  
  
  
  
    /*========================= Left foot ===============================*/
    var left_foot = generateEllipsoid(0, 0, 0, 0.06, 0.27, 0.06, 50, [0.4, 0.4, 0.4]);
    console.log(left_foot.vertices); // Check the generated vertices
  
  
    var LEFT_FOOT_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_FOOT_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(left_foot.vertices), GL.STATIC_DRAW);
  
  
    var LEFT_FOOT_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, LEFT_FOOT_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(left_foot.indices), GL.STATIC_DRAW);
  
  
    var LEFT_FOOT_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_FOOT_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(left_foot.colors), GL.STATIC_DRAW);
    /*====================================================================== */
  
  
  
  
  
  
    /*========================= Right foot ===============================*/
    var right_foot = generateEllipsoid(0, 0, 0, 0.06, 0.27, 0.06, 50, [0.4, 0.4, 0.4]);
    console.log(right_foot.vertices); // Check the generated vertices
  
  
    var RIGHT_FOOT_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_FOOT_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(right_foot.vertices), GL.STATIC_DRAW);
  
  
    var RIGHT_FOOT_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, RIGHT_FOOT_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(right_foot.indices), GL.STATIC_DRAW);
  
  
    var RIGHT_FOOT_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_FOOT_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(right_foot.colors), GL.STATIC_DRAW);
    /*====================================================================== */
  
  
  
  
  
  
    /*========================= Tail ===============================*/
    var tail = generateEllipsoid(0, 0, 0, 0.02, 0.3, 0.02, 50, [0.4, 0.4, 0.4]);
    console.log(tail.vertices); // Check the generated vertices
  
  
    var TAIL_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, TAIL_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tail.vertices), GL.STATIC_DRAW);
  
  
    var TAIL_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, TAIL_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(tail.indices), GL.STATIC_DRAW);
  
  
    var TAIL_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, TAIL_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tail.colors), GL.STATIC_DRAW);
    /*====================================================================== */
  
  
  
  
  
  
    /*========================= NOSE ===============================*/
    var nose = generateEllipsoid(0, 0, 0, 0.2, 0.3, 0.2, 50, [0.5, 0.5, 1.0]);
    console.log(nose.vertices); // Check the generated vertices
  
  
    var NOSE_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, NOSE_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(nose.vertices), GL.STATIC_DRAW);
  
  
    var NOSE_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, NOSE_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(nose.indices), GL.STATIC_DRAW);
  
  
    var NOSE_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, NOSE_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(nose.colors), GL.STATIC_DRAW);
  
  
    var nose_2 = generateEllipsoid(0, 0, 0, 0.07, 0.05, 0.05, 50, [1.0, 0.0, 0.7]);
    console.log(nose_2.vertices); // Check the generated vertices
  
  
    var NOSE_VERTEX_2 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, NOSE_VERTEX_2);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(nose_2.vertices), GL.STATIC_DRAW);
  
  
    var NOSE_INDICES_2 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, NOSE_INDICES_2);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(nose_2.indices), GL.STATIC_DRAW);
  
  
    var NOSE_COLOR_2 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, NOSE_COLOR_2);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(nose_2.colors), GL.STATIC_DRAW);
    /*====================================================================== */
  
  
  
  
  
  
    /*========================= Mouth ===============================*/
    var mainRadius = 0.3;
    var tubeRadius = 0.1;
    var numMainPoints = 60;
    var numTubePoints = 30;
    var torus_mulut_tv = halfbezierTorus(mainRadius, tubeRadius, numMainPoints, numTubePoints, 0.1, [0.0, 0.0, 0.0]);
  
  
    var MULUT_TV_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, MULUT_TV_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(torus_mulut_tv.vertices), GL.STATIC_DRAW);
  
  
    var MULUT_TV_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, MULUT_TV_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(torus_mulut_tv.indices), GL.STATIC_DRAW);
  
  
    var MULUT_TV_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, MULUT_TV_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(torus_mulut_tv.colors), GL.STATIC_DRAW);
    /*====================================================================== */
  
  
  
  
  
  
    /*========================= Left Ear ===============================*/
    var left_ear = generateEllipsoid(0, 0, 0, 0.06, 0.2, 0.1, 50, [0.6, 0.3, 1.0]);
    console.log(left_ear.vertices); // Check the generated vertices
  
  
    var LEFT_EAR_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_EAR_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(left_ear.vertices), GL.STATIC_DRAW);
  
  
    var LEFT_EAR_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, LEFT_EAR_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(left_ear.indices), GL.STATIC_DRAW);
  
  
    var LEFT_EAR_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_EAR_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(left_ear.colors), GL.STATIC_DRAW);
    /*====================================================================== */
  
  
  
  
  
  
    /*========================= Right Ear ===============================*/
    var right_ear = generateEllipsoid(0, 0, 0, 0.06, 0.2, 0.1, 50, [0.6, 0.3, 1.0]);
    console.log(right_ear.vertices); // Check the generated vertices
  
  
    var RIGHT_EAR_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_EAR_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(right_ear.vertices), GL.STATIC_DRAW);
  
  
    var RIGHT_EAR_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, RIGHT_EAR_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(right_ear.indices), GL.STATIC_DRAW);
  
  
    var RIGHT_EAR_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_EAR_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(right_ear.colors), GL.STATIC_DRAW);
    /*====================================================================== */
  
  
  
  
  
  
    /*========================= Janggut ========================== */
    var janggut_1 = createCylinder(GL, 0.01, 0.09, 50, 0, 0, 0);
    var janggut_2 = createCylinder(GL, 0.01, 0.09, 50, 0, 0, 0);
    var janggut_3 = createCylinder(GL, 0.01, 0.09, 50, 0, 0, 0);
    var janggut_4 = createCylinder(GL, 0.01, 0.09, 50, 0, 0, 0);
    var janggut_5 = createCylinder(GL, 0.01, 0.09, 50, 0, 0, 0);
    var janggut_6 = createCylinder(GL, 0.01, 0.09, 50, 0, 0, 0);
    var janggut_7 = createCylinder(GL, 0.01, 0.09, 50, 0, 0, 0);
  
  
    var JANGGUT_1_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_1_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(janggut_1.vertices), GL.STATIC_DRAW);
  
  
    var JANGGUT_1_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_1_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(janggut_1.colors), GL.STATIC_DRAW);
  
  
    var JANGGUT_2_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_2_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(janggut_2.vertices), GL.STATIC_DRAW);
  
  
    var JANGGUT_2_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_2_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(janggut_2.colors), GL.STATIC_DRAW);
  
  
    var JANGGUT_3_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_3_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(janggut_3.vertices), GL.STATIC_DRAW);
  
  
    var JANGGUT_3_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_3_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(janggut_3.colors), GL.STATIC_DRAW);
  
  
    var JANGGUT_4_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_4_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(janggut_4.vertices), GL.STATIC_DRAW);
  
  
    var JANGGUT_4_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_4_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(janggut_4.colors), GL.STATIC_DRAW);
  
  
    var JANGGUT_5_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_5_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(janggut_5.vertices), GL.STATIC_DRAW);
  
  
    var JANGGUT_5_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_5_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(janggut_5.colors), GL.STATIC_DRAW);
  
  
    var JANGGUT_6_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_6_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(janggut_6.vertices), GL.STATIC_DRAW);
  
  
    var JANGGUT_6_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_6_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(janggut_6.colors), GL.STATIC_DRAW);
  
  
    var JANGGUT_7_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_7_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(janggut_7.vertices), GL.STATIC_DRAW);
  
  
    var JANGGUT_7_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_7_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(janggut_7.colors), GL.STATIC_DRAW);
    /*====================================================================== */
  
  
  
  
  
  
    /*========================= Kalung_T.V ===============================*/
    var mainRadius = 5.4;
    var tubeRadius = 0.2;
    var numMainPoints = 60;
    var numTubePoints = 30;
    var torus_TV = bezierTorus(mainRadius, tubeRadius, numMainPoints, numTubePoints, 0.2, [1.0, 0.5, 0.0]);
  
  
    var TORUS_TV_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, TORUS_TV_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(torus_TV.vertices), GL.STATIC_DRAW);
  
  
    var TORUS_TV_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, TORUS_TV_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(torus_TV.indices), GL.STATIC_DRAW);
  
  
    var TORUS_TV_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, TORUS_TV_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(torus_TV.colors), GL.STATIC_DRAW);
    /*====================================================================== */
  
  
    //##############################################################################################################
    //##############################################################################################################
    //##############################################################################################################
  
  
  
  
  
  
    //##############################################################################################################
    //################################# GUNTER, Ice King's Minion ##################################################
    //##############################################################################################################
  
  
    /*========================= GUNTER BODY ===============================*/
    var gunter_body = generateElipticPara(1, 300, 1, 0, 0, 0, false);
  
  
    var GUNTER_BODY_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_BODY_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(gunter_body.vertices), GL.STATIC_DRAW);
  
  
    var GUNTER_BODY_FACES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, GUNTER_BODY_FACES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(gunter_body.indices), GL.STATIC_DRAW);
  
  
    var GUNTER_BODY_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_BODY_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(gunter_body.colors), GL.STATIC_DRAW);
    /*====================================================================== */
  
  
  
  
  
  
    /*========================= GUNTER BODY 2 ===============================*/
    var gunter_body2 = generateElipticPara(1, 300, 0.9, 1, 1, 1, false);
  
  
    var GUNTER_BODY_VERTEX2 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_BODY_VERTEX2);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(gunter_body2.vertices), GL.STATIC_DRAW);
  
  
    var GUNTER_BODY_FACES2 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, GUNTER_BODY_FACES2);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(gunter_body2.indices), GL.STATIC_DRAW);
  
  
    var GUNTER_BODY_COLOR2 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_BODY_COLOR2);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(gunter_body2.colors), GL.STATIC_DRAW);
    /*====================================================================== */
  
  
  
  
  
  
    /*========================= THE GUNTER_MOUTH ========================= */
    var GUNTER_MOUTH = generateCone(300, 0.3, 1.4, [1.0, 1.0, 0.0], 1); // Adjust parameters as needed
  
  
    var GUNTER_MOUTH_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_MOUTH_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(GUNTER_MOUTH.vertices), GL.STATIC_DRAW);
  
  
    var GUNTER_MOUTH_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_MOUTH_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(GUNTER_MOUTH.colors), GL.STATIC_DRAW);
    /*================================================================ */
  
  
  
  
  
  
    /*========================= THE Left GUNTER_EYE ========================= */
    var GUNTER_EYE = generateSphere(1, 300, 0.2);
    console.log(GUNTER_EYE.vertices);
  
  
    var GUNTER_EYE_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_EYE_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(GUNTER_EYE.vertices), GL.STATIC_DRAW);
  
  
    var GUNTER_EYE_FACES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, GUNTER_EYE_FACES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(GUNTER_EYE.indices), GL.STATIC_DRAW);
  
  
    var GUNTER_EYE_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_EYE_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(GUNTER_EYE.colors), GL.STATIC_DRAW);
    /*================================================================ */
  
  
  
  
  
  
    /*========================= The Right GUNTER_EYE ===============================*/
    var GUNTER_EYE_VERTEX2 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_EYE_VERTEX2);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(GUNTER_EYE.vertices), GL.STATIC_DRAW);
  
  
    var GUNTER_EYE_FACES2 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, GUNTER_EYE_FACES2);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(GUNTER_EYE.indices), GL.STATIC_DRAW);
  
  
    var GUNTER_EYE_COLOR2 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_EYE_COLOR2);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(GUNTER_EYE.colors), GL.STATIC_DRAW);
    /*================================================================ */
  
  
  
  
  
  
    /*========================= The Right GUNTER_HAND ===============================*/
    var GUNTER_HAND = generateEllipsoid(0, 0.7, 0.3, 0.15, 0.1, 0.6, 50, [0, 0, 0]);
  
  
    var GUNTER_HAND_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_HAND_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(GUNTER_HAND.vertices), GL.STATIC_DRAW);
  
  
    var GUNTER_HAND_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, GUNTER_HAND_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(GUNTER_HAND.indices), GL.STATIC_DRAW);
  
  
    var GUNTER_HAND_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_HAND_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(GUNTER_HAND.colors), GL.STATIC_DRAW);
    /*================================================================ */
  
  
  
  
  
  
    /*========================= The Left GUNTER_HAND ===============================*/
    var GUNTER_HAND_VERTEX2 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_HAND_VERTEX2);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(GUNTER_HAND.vertices), GL.STATIC_DRAW);
  
  
    var GUNTER_HAND_INDICES2 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, GUNTER_HAND_INDICES2);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(GUNTER_HAND.indices), GL.STATIC_DRAW);
  
  
    var GUNTER_HAND_COLOR2 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_HAND_COLOR2);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(GUNTER_HAND.colors), GL.STATIC_DRAW);
    /*================================================================ */
  
  
  
  
  
  
    /*========================= The Right GUNTER_FEET ===============================*/
    var GUNTER_FEET = generateHalfEllipsoid(0.2, 50, [1.0, 1.0, 0.0], [1.0, 1.0, 0.0], [0, 0, 0]);
  
  
    var GUNTER_FEET_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_FEET_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(GUNTER_FEET.vertices), GL.STATIC_DRAW);
  
  
    var GUNTER_FEET_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, GUNTER_FEET_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(GUNTER_FEET.indices), GL.STATIC_DRAW);
  
  
    var GUNTER_FEET_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_FEET_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(GUNTER_FEET.colors), GL.STATIC_DRAW);
    /*================================================================ */
  
  
  
  
  
  
    /*========================= The Right GUNTER_FEET ===============================*/
    var GUNTER_FEET_VERTEX2 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_FEET_VERTEX2);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(GUNTER_FEET.vertices), GL.STATIC_DRAW);
  
  
    var GUNTER_FEET_INDICES2 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, GUNTER_FEET_INDICES2);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(GUNTER_FEET.indices), GL.STATIC_DRAW);
  
  
    var GUNTER_FEET_COLOR2 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_FEET_COLOR2);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(GUNTER_FEET.colors), GL.STATIC_DRAW);
    /*================================================================ */
  
  
  
  
  
  
    /*========================= TORUS ===============================*/
    var mainRadius = 6.2;
    var tubeRadius = 0.2;
    var numMainPoints = 20;
    var numTubePoints = 10;
    var GUNTER_TORUS = bezierTorus(mainRadius, tubeRadius, numMainPoints, numTubePoints, 0.2, [1.0, 1.0, 0.0]);
  
  
    var GUNTER_TORUS_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_TORUS_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(GUNTER_TORUS.vertices), GL.STATIC_DRAW);
  
  
    var GUNTER_TORUS_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, GUNTER_TORUS_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(GUNTER_TORUS.indices), GL.STATIC_DRAW);
  
  
    var GUNTER_TORUS_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_TORUS_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(GUNTER_TORUS.colors), GL.STATIC_DRAW);
    //===================================================================
    
  
  
    /*=========================  permata gunter ===============================*/
    var permata = createCylinderJake(GL, 0.2, 0.1, 50, 1, 1, 0);
    var PERMATA_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, PERMATA_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(permata.vertices), GL.STATIC_DRAW);
  
  
    var PERMATA_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, PERMATA_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(permata.colors), GL.STATIC_DRAW);
  
  
    //##############################################################################################################
    //##############################################################################################################
    //##############################################################################################################
  
  
  
  
  
  
    //##############################################################################################################
    //############################### JAKE The Dog #################################################################
    //##############################################################################################################
  
  
    //=========================== BODY ====================================/
    //--------------------------- Middle ---------------------------------/
    var bodyCylinder_jake = createCylinderJake(GL, 1, 2, 50, 251 / 255, 215 / 255, 91 / 255);
    console.log(bodyCylinder_jake.vertices);
  
  
    var BODY_CYLINDER_VERTEX_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BODY_CYLINDER_VERTEX_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(bodyCylinder_jake.vertices), GL.STATIC_DRAW);
  
  
    var BODY_CYLINDER_COLOR_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BODY_CYLINDER_COLOR_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(bodyCylinder_jake.colors), GL.STATIC_DRAW);
    //----------------------------------------------------------------------/
    //---------------------------- Top -----------------------------------/
    var bodytop_jake = generateElipticPara(1.35, 350, 1, 251 / 255, 215 / 255, 91 / 255, false);
    console.log(bodytop_jake.vertices);
  
  
    var BODY_TOP_VERTEX_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BODY_TOP_VERTEX_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(bodytop_jake.vertices), GL.STATIC_DRAW);
  
  
    var BODY_TOP_INDICES_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, BODY_TOP_INDICES_JAKE);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(bodytop_jake.indices), GL.STATIC_DRAW);
  
  
    var BODY_TOP_COLOR_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BODY_TOP_COLOR_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(bodytop_jake.colors), GL.STATIC_DRAW);
    //----------------------------------------------------------------------/
    //---------------------------- Bottom --------------------------------/
    var bodybottom_jake = generateElipticPara(1.35, 350, 1, 251 / 255, 215 / 255, 91 / 255, false);
    console.log(bodybottom_jake.vertices);
  
  
    var BODY_BOTTOM_VERTEX_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BODY_BOTTOM_VERTEX_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(bodybottom_jake.vertices), GL.STATIC_DRAW);
  
  
    var BODY_BOTTOM_INDICES_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, BODY_BOTTOM_INDICES_JAKE);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(bodybottom_jake.indices), GL.STATIC_DRAW);
  
  
    var BODY_BOTTOM_COLOR_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BODY_BOTTOM_COLOR_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(bodybottom_jake.colors), GL.STATIC_DRAW);
    //----------------------------------------------------------------------/
    //==========================================================================/
  
  
  
  
  
  
    //========================== ARM ====================================/
    //-------------------------- Left -------------------------------------/
    var numTri = 50;
    var leftarmcylinder_jake = createCylinderJake(GL, 0.125, 2, numTri, 251 / 255, 215 / 255, 91 / 255);
    console.log(leftarmcylinder_jake.vertices);
  
  
    var LEFT_ARM_CYLINDER_VERTEX_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_ARM_CYLINDER_VERTEX_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(leftarmcylinder_jake.vertices), GL.STATIC_DRAW);
  
  
    var LEFT_ARM_CYLINDER_COLOR_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_ARM_CYLINDER_COLOR_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(leftarmcylinder_jake.colors), GL.STATIC_DRAW);
    //----------------------------------------------------------------------/
    //-------------------------- Right ------------------------------------/
    var rightarmcylinder_jake = createCylinderJake(GL, 0.125, 2, 50, 251 / 255, 215 / 255, 91 / 255);
    console.log(rightarmcylinder_jake.vertices);
  
  
    var RIGHT_ARM_CYLINDER_VERTEX_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_ARM_CYLINDER_VERTEX_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rightarmcylinder_jake.vertices), GL.STATIC_DRAW);
  
  
    var RIGHT_ARM_CYLINDER_COLOR_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_ARM_CYLINDER_COLOR_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rightarmcylinder_jake.colors), GL.STATIC_DRAW);
    //----------------------------------------------------------------------/
    //=====================================================================/
  
  
  
  
  
  
    //========================== HAND ====================================/
    //--------------------------- Left -----------------------------------/
    var lefthand_jake = generateElipticPara(0.14, 50, 1, 251 / 255, 215 / 255, 91 / 255, false);
    console.log(lefthand_jake.vertices);
  
  
    var LEFT_HAND_VERTEX_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_HAND_VERTEX_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(lefthand_jake.vertices), GL.STATIC_DRAW);
  
  
    var LEFT_HAND_INDICES_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, LEFT_HAND_INDICES_JAKE);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(lefthand_jake.indices), GL.STATIC_DRAW);
  
  
    var LEFT_HAND_COLOR_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_HAND_COLOR_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(lefthand_jake.colors), GL.STATIC_DRAW);
    //----------------------------------------------------------------------/
    //---------------------------- Right -----------------------------------/
    var righthand_jake = generateElipticPara(0.14, 50, 1, 251 / 255, 215 / 255, 91 / 255, false);
    console.log(righthand_jake.vertices);
  
  
    var RIGHT_HAND_VERTEX_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_HAND_VERTEX_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(righthand_jake.vertices), GL.STATIC_DRAW);
  
  
    var RIGHT_HAND_INDICES_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, RIGHT_HAND_INDICES_JAKE);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(righthand_jake.indices), GL.STATIC_DRAW);
  
  
    var RIGHT_HAND_COLOR_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_HAND_COLOR_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(righthand_jake.colors), GL.STATIC_DRAW);
    //----------------------------------------------------------------------/
    //====================================================================/
  
  
  
  
  
  
    //========================== FEET ====================================/
    //-------------------------- Left ------------------------------------/
    var leftfeetcylinder_jake = createCylinder(GL, 0.125, 2.5, 50, 251 / 255, 215 / 255, 91 / 255);
    console.log(leftfeetcylinder_jake.vertices);
  
  
    var LEFT_FEET_CYLINDER_VERTEX_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_FEET_CYLINDER_VERTEX_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(leftfeetcylinder_jake.vertices), GL.STATIC_DRAW);
  
  
    var LEFT_FEET_CYLINDER_COLOR_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_FEET_CYLINDER_COLOR_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(leftfeetcylinder_jake.colors), GL.STATIC_DRAW);
    //----------------------------------------------------------------------/
    //-------------------------- Right ------------------------------------/
    var rightfeetcylinder_jake = createCylinder(GL, 0.125, 2.5, 50, 251 / 255, 215 / 255, 91 / 255);
    console.log(rightfeetcylinder_jake.vertices);
  
  
    var RIGHT_FEET_CYLINDER_VERTEX_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_FEET_CYLINDER_VERTEX_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rightfeetcylinder_jake.vertices), GL.STATIC_DRAW);
  
  
    var RIGHT_FEET_CYLINDER_COLOR_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_FEET_CYLINDER_COLOR_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rightfeetcylinder_jake.colors), GL.STATIC_DRAW);
    //----------------------------------------------------------------------/
    //====================================================================/
  
  
  
  
  
  
    //========================== FOOT ====================================/
    //---------------------------- Left -----------------------------------/
    var leftfoot_jake = generateElipticPara(0.14, 50, 1, 251 / 255, 215 / 255, 91 / 255, false);
    console.log(leftfoot_jake.vertices);
  
  
    var LEFT_FOOT_VERTEX_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_FOOT_VERTEX_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(leftfoot_jake.vertices), GL.STATIC_DRAW);
  
  
    var LEFT_FOOT_INDICES_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, LEFT_FOOT_INDICES_JAKE);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(leftfoot_jake.indices), GL.STATIC_DRAW);
  
  
    var LEFT_FOOT_COLOR_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_FOOT_COLOR_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(leftfoot_jake.colors), GL.STATIC_DRAW);
    //----------------------------------------------------------------------/
    //---------------------------- Right -----------------------------------/
    var rightfoot_jake = generateElipticPara(0.14, 50, 1, 251 / 255, 215 / 255, 91 / 255, false);
    console.log(rightfoot_jake.vertices);
  
  
    var RIGHT_FOOT_VERTEX_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_FOOT_VERTEX_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rightfoot_jake.vertices), GL.STATIC_DRAW);
  
  
    var RIGHT_FOOT_INDICES_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, RIGHT_FOOT_INDICES_JAKE);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rightfoot_jake.indices), GL.STATIC_DRAW);
  
  
    var RIGHT_FOOT_COLOR_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_FOOT_COLOR_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rightfoot_jake.colors), GL.STATIC_DRAW);
    //----------------------------------------------------------------------/
    //====================================================================/
  
  
  
  
  
  
    //============================= EYE ====================================/
    //-------------------------- Left -------------------------------------/
    var lefteye_jake = createCylinder(GL, 0.3, -0.1, 50, 0, 0, 0);
    console.log(lefteye_jake.vertices);
  
  
    var LEFT_EYE_VERTEX_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_EYE_VERTEX_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(lefteye_jake.vertices), GL.STATIC_DRAW);
  
  
    var LEFT_EYE_COLOR_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_EYE_COLOR_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(lefteye_jake.colors), GL.STATIC_DRAW);
    //----------------------------------------------------------------------/
    //-------------------------- Right ------------------------------------/
    var righteye_jake = createCylinder(GL, 0.3, -0.1, 50, 0, 0, 0);
    console.log(righteye_jake.vertices);
  
  
    var RIGHT_EYE_VERTEX_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_EYE_VERTEX_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(righteye_jake.vertices), GL.STATIC_DRAW);
  
  
    var RIGHT_EYE_COLOR_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_EYE_COLOR_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(righteye_jake.colors), GL.STATIC_DRAW);
    //----------------------------------------------------------------------/
    //=======================================================================/
  
  
  
  
  
  
    //============================= NOSE ====================================/
    //-------------------------- BLACK ONE ----------------------------------/
    var blacknose_jake = generateEllipsoid(0, 0, 0, 0.07, 0.05, 0.05, 50, [0.0, 0.0, 0.0]);
  
  
    var BLACK_NOSE_VERTEX_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BLACK_NOSE_VERTEX_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(blacknose_jake.vertices), GL.STATIC_DRAW);
  
  
    var BLACK_NOSE_INDICES_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, BLACK_NOSE_INDICES_JAKE);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(blacknose_jake.indices), GL.STATIC_DRAW);
  
  
    var BLACK_NOSE_COLOR_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BLACK_NOSE_COLOR_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(blacknose_jake.colors), GL.STATIC_DRAW);
    //----------------------------------------------------------------------/
    //--------------------------- CURVE TORUS ------------------------------/
    var unose_torus_jake = halfbezierTorus(2, 1, 20, 10, 0.1, [255 / 255, 255 / 255, 0 / 255]);
  
  
    var UNOSE_TORUS_VERTEX_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, UNOSE_TORUS_VERTEX_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(unose_torus_jake.vertices), GL.STATIC_DRAW);
  
  
    var UNOSE_TORUS_INDICES_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, UNOSE_TORUS_INDICES_JAKE);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(unose_torus_jake.indices), GL.STATIC_DRAW);
  
  
    var UNOSE_TORUS_COLOR_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, UNOSE_TORUS_COLOR_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(unose_torus_jake.colors), GL.STATIC_DRAW);
    //----------------------------------------------------------------------/
    //=======================================================================/
  
  
  
  
  
  
    //============================== EAR =====================================/
    //--------------------------- Left -----------------------------------/
    var leftear_jake = generateElipticPara(0.14, 50, 1, 251 / 255, 215 / 255, 91 / 255, false);
    console.log(leftear_jake.vertices);
  
  
    var LEFT_EAR_VERTEX_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_EAR_VERTEX_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(leftear_jake.vertices), GL.STATIC_DRAW);
  
  
    var LEFT_EAR_INDICES_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, LEFT_EAR_INDICES_JAKE);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(leftear_jake.indices), GL.STATIC_DRAW);
  
  
    var LEFT_EAR_COLOR_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_EAR_COLOR_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(leftear_jake.colors), GL.STATIC_DRAW);
    //----------------------------------------------------------------------/
    //---------------------------- Right -----------------------------------/
    var rightear_jake = generateElipticPara(0.14, 50, 1, 251 / 255, 215 / 255, 91 / 255, false);
    console.log(rightear_jake.vertices);
  
  
    var RIGHT_EAR_VERTEX_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_EAR_VERTEX_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rightear_jake.vertices), GL.STATIC_DRAW);
  
  
    var RIGHT_EAR_INDICES_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, RIGHT_EAR_INDICES_JAKE);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(rightear_jake.indices), GL.STATIC_DRAW);
  
  
    var RIGHT_EAR_COLOR_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_EAR_COLOR_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(rightear_jake.colors), GL.STATIC_DRAW);
    //----------------------------------------------------------------------/
    //========================================================================/
  
  
  
  
  
  
    //========================== TAIL ====================================/
    var tail_jake = generateElipticPara(0.14, 50, 1, 251 / 255, 215 / 255, 91 / 255, false);
    console.log(tail_jake.vertices);
  
  
    var TAIL_VERTEX_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, TAIL_VERTEX_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tail_jake.vertices), GL.STATIC_DRAW);
  
  
    var TAIL_INDICES_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, TAIL_INDICES_JAKE);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(tail_jake.indices), GL.STATIC_DRAW);
  
  
    var TAIL_COLOR_JAKE = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, TAIL_COLOR_JAKE);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(tail_jake.colors), GL.STATIC_DRAW);
    //========================================================================/
  
  
    //##############################################################################################################
    //##############################################################################################################
    //##############################################################################################################
  
  
  
  
  
  
    //##############################################################################################################
    //########################## ENVIRONTMENT OF ADVENTURE TIME ####################################################
    //##############################################################################################################
  
  
    //================================== BASE  ====================================/
    var base = generateHalfEllipsoid(30, 70, [0.3, 0.3, 0.3], [1, 1, 1], [0, 0.9, 0]);
  
  
    var BASE_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BASE_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(base.vertices), GL.STATIC_DRAW);
  
  
    var BASE_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, BASE_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(base.indices), GL.STATIC_DRAW);
  
  
    var BASE_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BASE_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(base.colors), GL.STATIC_DRAW);
    //==============================================================================/
  
  
  
  
  
  
    //================================== BASE 2 ====================================/
    var base2 = generateHalfEllipsoid(30, 70, [0.3, 0.3, 0.3], [1, 1, 1], [0, 0, 1]);
  
  
    var BASE2_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BASE2_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(base2.vertices), GL.STATIC_DRAW);
  
  
    var BASE2_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, BASE2_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(base2.indices), GL.STATIC_DRAW);
  
  
    var BASE2_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BASE2_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(base2.colors), GL.STATIC_DRAW);
    //==============================================================================/
  
  
  
  
  
  
     //================================== BASE 3 ====================================/
     var base3 = generateHalfEllipsoid(13, 70, [0.3, 0.3, 0.3], [1, 1, 1], [0.6, 0.3, 0.0]);
  
  
     var BASE3_VERTEX = GL.createBuffer();
     GL.bindBuffer(GL.ARRAY_BUFFER, BASE3_VERTEX);
     GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(base3.vertices), GL.STATIC_DRAW);
   
     var BASE3_INDICES = GL.createBuffer();
     GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, BASE3_INDICES);
     GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(base3.indices), GL.STATIC_DRAW);
   
     var BASE3_COLOR = GL.createBuffer();
     GL.bindBuffer(GL.ARRAY_BUFFER, BASE3_COLOR);
     GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(base3.colors), GL.STATIC_DRAW);
     //==============================================================================/
  
  
  
  
  
  
    //================================= BATANG MAIN TREE HOUSE =========================================/
    var batang = generateCone(300, 0.1, 12, [0.5, 0.2, 0], 70);
  
  
    var BATANG_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BATANG_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(batang.vertices), GL.STATIC_DRAW);
  
  
    var BATANG_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, BATANG_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(batang.indices), GL.STATIC_DRAW);
  
  
    var BATANG_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BATANG_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(batang.colors), GL.STATIC_DRAW);
  
  
    var batang1 = generateCone(300, 0.1, 7, [0.4, 0.13, 0.05], 10);
  
  
    var BATANG1_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BATANG1_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(batang1.vertices), GL.STATIC_DRAW);
  
  
    var BATANG1_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, BATANG1_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(batang1.indices), GL.STATIC_DRAW);
  
  
    var BATANG1_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BATANG1_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(batang1.colors), GL.STATIC_DRAW);
  
  
    var batang2 = generateCone(300, 0.2, 10, [0.4, 0.13, 0.05], 10);
  
  
    var BATANG2_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BATANG2_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(batang2.vertices), GL.STATIC_DRAW);
  
  
    var BATANG2_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, BATANG2_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(batang2.indices), GL.STATIC_DRAW);
  
  
    var BATANG2_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BATANG2_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(batang2.colors), GL.STATIC_DRAW);
  
  
    var batang3 = generateCone(300, 0.01, 9, [0.5, 0.2, 0], 70);
  
  
    var BATANG3_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BATANG3_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(batang3.vertices), GL.STATIC_DRAW);
  
  
    var BATANG3_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, BATANG3_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(batang3.indices), GL.STATIC_DRAW);
  
  
    var BATANG3_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, BATANG3_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(batang3.colors), GL.STATIC_DRAW);
    //==============================================================================================/
  
  
  
  
  
  
    //========================================= DAUN BESAR ==========================================/
    //------------------------- Daun Besar ------------------------------------/
    var daun = generateElipticPara(1, 100, 7, 0, 0.5, 0, true);
  
  
    var DAUN_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, DAUN_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(daun.vertices), GL.STATIC_DRAW);
  
  
    var DAUN_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, DAUN_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(daun.indices), GL.STATIC_DRAW);
  
  
    var DAUN_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, DAUN_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(daun.colors), GL.STATIC_DRAW);
    //--------------------------------------------------------------------------/
  
  
    //------------------------- Daun Besar 1 ------------------------------------/
    var daun1 = generateElipticPara(1, 100, 3, 0, 0.3, 0, true);
  
  
    var DAUN1_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, DAUN1_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(daun1.vertices), GL.STATIC_DRAW);
  
  
    var DAUN1_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, DAUN1_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(daun1.indices), GL.STATIC_DRAW);
  
  
    var DAUN1_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, DAUN1_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(daun1.colors), GL.STATIC_DRAW);
    //----------------------------------------------------------------------------/
  
  
    //------------------------- Daun Besar 2 ------------------------------------/
    var daun2 = generateElipticPara(1, 100, 3, 0, 0.7, 0, true);
  
  
    var DAUN2_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, DAUN2_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(daun2.vertices), GL.STATIC_DRAW);
  
  
    var DAUN2_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, DAUN2_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(daun2.indices), GL.STATIC_DRAW);
  
  
    var DAUN2_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, DAUN2_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(daun2.colors), GL.STATIC_DRAW);
    //-----------------------------------------------------------------------------/
  
  
    //---------------------------- Daun Besar 3 ------------------------------------/
    var daun3 = generateElipticPara(1, 100, 3.5, 0, 0.6, 0, true);
  
  
    var DAUN3_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, DAUN3_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(daun3.vertices), GL.STATIC_DRAW);
  
  
    var DAUN3_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, DAUN3_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(daun3.indices), GL.STATIC_DRAW);
  
  
    var DAUN3_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, DAUN3_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(daun3.colors), GL.STATIC_DRAW);
    //------------------------------------------------------------------------------/
  //---------------------------------------------ICE KINGDOM------------------------------/
    const numSegments = 10;
    const width = 6.0;
    const height = 8.0;
    const numZigs = 10;
    
    const zig = generateConnectedZigzagPath(numSegments, width, height, numZigs, 5);
  
  
  
  
    var ZIG_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, ZIG_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(zig.vertices), GL.STATIC_DRAW);
  
  
    var ZIG_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, ZIG_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(zig.indices), GL.STATIC_DRAW);
  
  
    var ZIG_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, ZIG_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(zig.colors), GL.STATIC_DRAW);
  
  
    var mount = generateMount(20, 16, 60,  [1.0, 1.0, 1.0], [0.1, 60/255 , 249/255], 1); 
  
  
    var MOUNT_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER,MOUNT_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(mount.vertices), GL.STATIC_DRAW);
  
  
    var MOUNT_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, MOUNT_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(mount.colors), GL.STATIC_DRAW);
  
  
    var mount1 = generateMount(20, 13, 40,  [1.0, 1.0, 1.0], [0.1, 140/255 , 249/255], 1); 
  
  
    var MOUNT_VERTEX1 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER,MOUNT_VERTEX1);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(mount1.vertices), GL.STATIC_DRAW);
  
  
    var MOUNT_COLOR1 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, MOUNT_COLOR1);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(mount1.colors), GL.STATIC_DRAW);
  
  
    var gua = generateMount(20, 2, 5,  [0.4, 0.4, 0.4], [0, 0, 0], 1); 
  
  
    var GUA_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER,GUA_VERTEX);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(gua.vertices), GL.STATIC_DRAW);
  
  
    var GUA_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, GUA_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(gua.colors), GL.STATIC_DRAW);
  
  
    var eye_kingdom = generateEllipsoid(0, 0, 0, 1, 1, 1, 50, [0, 0, 0]);
    var EYE_KINGDOM_VERTEX = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER,EYE_KINGDOM_VERTEX );
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(eye_kingdom.vertices), GL.STATIC_DRAW);
  
  
    var EYE_KINGDOM_INDICES = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, EYE_KINGDOM_INDICES);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(eye_kingdom.indices), GL.STATIC_DRAW);
  
  
    var  EYE_KINGDOM_COLOR = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, EYE_KINGDOM_COLOR);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(eye_kingdom.colors), GL.STATIC_DRAW);
    
    var EYE_KINGDOM_VERTEX1 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER,EYE_KINGDOM_VERTEX1 );
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(eye_kingdom.vertices), GL.STATIC_DRAW);
  
  
    var EYE_KINGDOM_INDICES1 = GL.createBuffer();
    GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, EYE_KINGDOM_INDICES1);
    GL.bufferData(GL.ELEMENT_ARRAY_BUFFER, new Uint16Array(eye_kingdom.indices), GL.STATIC_DRAW);
  
  
    var  EYE_KINGDOM_COLOR1 = GL.createBuffer();
    GL.bindBuffer(GL.ARRAY_BUFFER, EYE_KINGDOM_COLOR1);
    GL.bufferData(GL.ARRAY_BUFFER, new Float32Array(eye_kingdom.colors), GL.STATIC_DRAW);
  
  
    //==============================================================================================/
    //##############################################################################################################
    //##############################################################################################################
    //##############################################################################################################
  
  
  
  
  
  
    //##############################################################################################################
    //##############################################################################################################
    //##############################################################################################################
  
  
  
  
  
  
    //##############################################################################################################
    //##############################################################################################################
    //##############################################################################################################
  
  
  
  
  
  
    //##############################################################################################################
    //############################ MATRIX OF EACH CHARACTER ########################################################
    //##############################################################################################################
  
  
    /*========================= MATRIX T.V The Son of Jake ===============================*/
    var PROJECTION_MATRIX = LIBS.get_projection(40, CANVAS.width / CANVAS.height, 1, 100);
    var VIEW_MATRIX = LIBS.get_I4();
    var MODEL_MATRIX = LIBS.get_I4();
  
  
    var CYLINDER_MODEL_MATRIX = LIBS.get_I4();
    var CYLINDER_MODEL_MATRIX_2 = LIBS.get_I4();
  
  
    var Right_BROW_MODEL_MATRIX = LIBS.get_I4();
    var Left_BROW_MODEL_MATRIX = LIBS.get_I4();
  
  
    var Left_HAND_MODEL_MATRIX = LIBS.get_I4();
    var Right_HAND_MODEL_MATRIX = LIBS.get_I4();
  
  
    var Left_FOOT_MODEL_MATRIX = LIBS.get_I4();
    var Right_FOOT_MODEL_MATRIX = LIBS.get_I4();
  
  
    var TAIL_MODEL_MATRIX = LIBS.get_I4();
    var NOSE_MODEL_MATRIX = LIBS.get_I4();
    var NOSE_MODEL_MATRIX_2 = LIBS.get_I4();
    var MULUT_MODEL_MATRIX = LIBS.get_I4();
  
  
    var LEFT_EAR_MODEL_MATRIX = LIBS.get_I4();
    var RIGHT_EAR_MODEL_MATRIX = LIBS.get_I4();
  
  
    var JANGGUT_1_MODEL_MATRIX = LIBS.get_I4();
    var JANGGUT_2_MODEL_MATRIX = LIBS.get_I4();
    var JANGGUT_3_MODEL_MATRIX = LIBS.get_I4();
    var JANGGUT_4_MODEL_MATRIX = LIBS.get_I4();
    var JANGGUT_5_MODEL_MATRIX = LIBS.get_I4();
    var JANGGUT_6_MODEL_MATRIX = LIBS.get_I4();
    var JANGGUT_7_MODEL_MATRIX = LIBS.get_I4();
  
  
    var TORUS_TV_MODEL_MATRIX = LIBS.get_I4();
    /*====================================================================== */
  
  
  
  
  
  
    /*========================= MATRIX GUNTER ===============================*/
    var GUNTER_MODEL_MATRIX = LIBS.get_I4();
    var GUNTER_MODEL_MATRIX2 = LIBS.get_I4();
    var GUNTER_MOUTH_MODEL_MATRIX = LIBS.get_I4();
    var GUNTER_EYE_MODEL_MATRIX = LIBS.get_I4();
    var GUNTER_EYE_MODEL_MATRIX2 = LIBS.get_I4();
    var GUNTER_HAND_MODEL_MATRIX = LIBS.get_I4();
    var GUNTER_HAND_MODEL_MATRIX2 = LIBS.get_I4();
    var GUNTER_FEET_MODEL_MATRIX = LIBS.get_I4();
    var GUNTER_FEET_MODEL_MATRIX2 = LIBS.get_I4();
    var GUNTER_TORUS_MODEL_MATRIX = LIBS.get_I4();
    var GUNTER_PERMATA_MODEL_MATRIX = LIBS.get_I4();
    /*====================================================================== */
  
  
  
  
  
  
    //=========================== MATRIX JAKE ====================================/
    var PROJECTION_MATRIX_JAKE = LIBS.get_projection(40, CANVAS.width / CANVAS.height, 1, 100);
    var VIEW_MATRIX_JAKE = LIBS.get_I4();
    var MODEL_MATRIX_JAKE = LIBS.get_I4();
  
  
    var MODEL_BODY_CYLINDER_MATRIX_JAKE = LIBS.get_I4();
    var MODEL_BODY_TOP_MATRIX_JAKE = LIBS.get_I4();
    var MODEL_BODY_BOTTOM_MATRIX_JAKE = LIBS.get_I4();
  
  
    var MODEL_LEFT_ARM_MATRIX_JAKE = LIBS.get_I4();
    var MODEL_RIGHT_ARM_MATRIX_JAKE = LIBS.get_I4();
  
  
    var MODEL_LEFT_FEET_MATRIX_JAKE = LIBS.get_I4();
    var MODEL_RIGHT_FEET_MATRIX_JAKE = LIBS.get_I4();
  
  
    var MODEL_LEFT_FOOT_MATRIX_JAKE = LIBS.get_I4();
    var MODEL_RIGHT_FOOT_MATRIX_JAKE = LIBS.get_I4();
  
  
    var MODEL_LEFT_EYE_MATRIX_JAKE = LIBS.get_I4();
    var MODEL_RIGHT_EYE_MATRIX_JAKE = LIBS.get_I4();
  
  
    var MODEL_BLACK_NOSE_MATRIX_JAKE = LIBS.get_I4();
    var MODEL_CURVE_TORUS_MATRIX_JAKE = LIBS.get_I4();
  
  
    var MODEL_LEFT_EAR_MATRIX_JAKE = LIBS.get_I4();
    var MODEL_RIGHT_EAR_MATRIX_JAKE = LIBS.get_I4();
  
  
    var MODEL_TAIL_MATRIX_JAKE = LIBS.get_I4();
    //=================================================================================/
  
  
  
  
    //=========================== MATRIX ENVIRONMENT ====================================/
    var BASE_MODEL_MATRIX = LIBS.get_I4();
    var BASE2_MODEL_MATRIX = LIBS.get_I4();
    var BASE3_MODEL_MATRIX = LIBS.get_I4();
  
  
    var BATANG_MODEL_MATRIX = LIBS.get_I4();
    var BATANG1_MODEL_MATRIX = LIBS.get_I4();
    var BATANG2_MODEL_MATRIX = LIBS.get_I4();
    var BATANG3_MODEL_MATRIX = LIBS.get_I4();
  
  
    var DAUN_MODEL_MATRIX = LIBS.get_I4();
    var DAUN1_MODEL_MATRIX = LIBS.get_I4();
    var DAUN2_MODEL_MATRIX = LIBS.get_I4();
    var DAUN3_MODEL_MATRIX = LIBS.get_I4();
  
  
    var ZIG_MODEL_MATRIX = LIBS.get_I4();
    var MOUNT_MODEL_MATRIX = LIBS.get_I4();
    var MOUNT_MODEL_MATRIX1 = LIBS.get_I4();
    var GUA_MODEL_MATRIX = LIBS.get_I4();
    var EYE_KINGDOM_MODEL_MATRIX = LIBS.get_I4();
    var EYE_KINGDOM_MODEL_MATRIX1= LIBS.get_I4();
    //=================================================================================/
  
  
  
  
  
  
  
  
  
  
  
  
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SETTING T.V $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    /*========================= SETTING THE LEFT EYE ========================== */
    LIBS.translateX(CYLINDER_MODEL_MATRIX, LIBS.degToRad(20)); // Move the cylinder along the X-axis
    LIBS.translateY(CYLINDER_MODEL_MATRIX, LIBS.degToRad(20)); // Move the cylinder along the Y-axis
    LIBS.translateZ(CYLINDER_MODEL_MATRIX, LIBS.degToRad(53)); // Move the cylinder along the Z-axis
    LIBS.rotateX(CYLINDER_MODEL_MATRIX, LIBS.degToRad(90)); // Rotate the cylinder around the X-axis
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING THE RIGHT EYE ========================== */
    LIBS.translateX(CYLINDER_MODEL_MATRIX_2, LIBS.degToRad(-20)); // Move the cylinder along the X-axis
    LIBS.translateY(CYLINDER_MODEL_MATRIX_2, LIBS.degToRad(20)); // Move the cylinder along the Y-axis
    LIBS.translateZ(CYLINDER_MODEL_MATRIX_2, LIBS.degToRad(53)); // Move the cylinder along the Z-axis
    LIBS.rotateX(CYLINDER_MODEL_MATRIX_2, LIBS.degToRad(90)); // Rotate the cylinder around the X-axis
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING THE RIGHT BROW ========================== */
    LIBS.translateX(Right_BROW_MODEL_MATRIX, LIBS.degToRad(-20)); // Move the cylinder along the X-axis
    LIBS.translateY(Right_BROW_MODEL_MATRIX, LIBS.degToRad(30)); // Move the cylinder along the Y-axis
    LIBS.translateZ(Right_BROW_MODEL_MATRIX, LIBS.degToRad(46)); // Move the cylinder along the Z-axis
    LIBS.rotateX(Right_BROW_MODEL_MATRIX, LIBS.degToRad(90)); // Rotate the cylinder around the X-axis
    LIBS.rotateZ(Right_BROW_MODEL_MATRIX, LIBS.degToRad(-9)); // Rotate the cylinder around the Y-axis
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING THE LEFT BROW ========================== */
    LIBS.translateX(Left_BROW_MODEL_MATRIX, LIBS.degToRad(20)); // Move the cylinder along the X-axis
    LIBS.translateY(Left_BROW_MODEL_MATRIX, LIBS.degToRad(30)); // Move the cylinder along the Y-axis
    LIBS.translateZ(Left_BROW_MODEL_MATRIX, LIBS.degToRad(46)); // Move the cylinder along the Z-axis
    LIBS.rotateX(Left_BROW_MODEL_MATRIX, LIBS.degToRad(90)); // Rotate the cylinder around the X-axis
    LIBS.rotateZ(Left_BROW_MODEL_MATRIX, LIBS.degToRad(10)); // Rotate the cylinder around the Y-axis
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING THE LEFT HAND ========================== */
    LIBS.translateX(Left_HAND_MODEL_MATRIX, LIBS.degToRad(47));
    LIBS.translateY(Left_HAND_MODEL_MATRIX, LIBS.degToRad(27));
    LIBS.rotateZ(Left_HAND_MODEL_MATRIX, LIBS.degToRad(20));
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING THE RIGHT HAND ========================== */
    LIBS.translateX(Right_HAND_MODEL_MATRIX, LIBS.degToRad(-47));
    LIBS.translateY(Right_HAND_MODEL_MATRIX, LIBS.degToRad(27));
    LIBS.rotateZ(Right_HAND_MODEL_MATRIX, LIBS.degToRad(-20));
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING THE LEFT FOOT ========================== */
    LIBS.translateX(Left_FOOT_MODEL_MATRIX, LIBS.degToRad(-20));
    LIBS.translateY(Left_FOOT_MODEL_MATRIX, LIBS.degToRad(-23));
    LIBS.rotateZ(Left_FOOT_MODEL_MATRIX, LIBS.degToRad(0));
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING THE RIGHT FOOT ========================== */
    LIBS.translateX(Right_FOOT_MODEL_MATRIX, LIBS.degToRad(20));
    LIBS.translateY(Right_FOOT_MODEL_MATRIX, LIBS.degToRad(-23));
    LIBS.rotateZ(Right_FOOT_MODEL_MATRIX, LIBS.degToRad(0));
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING THE TAIL ========================== */
    LIBS.translateX(TAIL_MODEL_MATRIX, LIBS.degToRad(0));
    LIBS.translateY(TAIL_MODEL_MATRIX, LIBS.degToRad(-10));
    LIBS.translateZ(TAIL_MODEL_MATRIX, LIBS.degToRad(-25));
    LIBS.rotateX(TAIL_MODEL_MATRIX, LIBS.degToRad(30));
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING THE NOSE ========================== */
    LIBS.translateX(NOSE_MODEL_MATRIX, LIBS.degToRad(0));
    LIBS.translateY(NOSE_MODEL_MATRIX, LIBS.degToRad(8));
    LIBS.translateZ(NOSE_MODEL_MATRIX, LIBS.degToRad(95));
    LIBS.rotateX(NOSE_MODEL_MATRIX, LIBS.degToRad(90));
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING THE NOSE 2 ========================== */
    LIBS.translateX(NOSE_MODEL_MATRIX_2, LIBS.degToRad(0));
    LIBS.translateY(NOSE_MODEL_MATRIX_2, LIBS.degToRad(11));
    LIBS.translateZ(NOSE_MODEL_MATRIX_2, LIBS.degToRad(110));
    LIBS.rotateX(NOSE_MODEL_MATRIX_2, LIBS.degToRad(90));
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING THE MOUTH ========================== */
    LIBS.translateX(MULUT_MODEL_MATRIX, LIBS.degToRad(0));
    LIBS.translateY(MULUT_MODEL_MATRIX, LIBS.degToRad(5));
    LIBS.translateZ(MULUT_MODEL_MATRIX, LIBS.degToRad(71.97));
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING THE LEFT EAR ========================== */
    LIBS.translateX(LEFT_EAR_MODEL_MATRIX, LIBS.degToRad(23));
    LIBS.translateY(LEFT_EAR_MODEL_MATRIX, LIBS.degToRad(45));
    LIBS.translateZ(LEFT_EAR_MODEL_MATRIX, LIBS.degToRad(9));
    LIBS.rotateZ(LEFT_EAR_MODEL_MATRIX, LIBS.degToRad(50));
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING THE RIGHT EAR ========================== */
    LIBS.translateX(RIGHT_EAR_MODEL_MATRIX, LIBS.degToRad(-23));
    LIBS.translateY(RIGHT_EAR_MODEL_MATRIX, LIBS.degToRad(45));
    LIBS.translateZ(RIGHT_EAR_MODEL_MATRIX, LIBS.degToRad(9));
    LIBS.rotateZ(RIGHT_EAR_MODEL_MATRIX, LIBS.degToRad(-50));
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING THE JANGGUT ========================== */
    LIBS.translateX(JANGGUT_1_MODEL_MATRIX, LIBS.degToRad(-20));
    LIBS.translateY(JANGGUT_1_MODEL_MATRIX, LIBS.degToRad(-7));
    LIBS.translateZ(JANGGUT_1_MODEL_MATRIX, LIBS.degToRad(55));
    LIBS.rotateX(JANGGUT_1_MODEL_MATRIX, LIBS.degToRad(-20));
  
  
    LIBS.translateX(JANGGUT_2_MODEL_MATRIX, LIBS.degToRad(-17));
    LIBS.translateY(JANGGUT_2_MODEL_MATRIX, LIBS.degToRad(1));
    LIBS.translateZ(JANGGUT_2_MODEL_MATRIX, LIBS.degToRad(56));
    LIBS.rotateX(JANGGUT_2_MODEL_MATRIX, LIBS.degToRad(-20));
  
  
    LIBS.translateX(JANGGUT_3_MODEL_MATRIX, LIBS.degToRad(16));
    LIBS.translateY(JANGGUT_3_MODEL_MATRIX, LIBS.degToRad(-2));
    LIBS.translateZ(JANGGUT_3_MODEL_MATRIX, LIBS.degToRad(56));
    LIBS.rotateX(JANGGUT_3_MODEL_MATRIX, LIBS.degToRad(-20));
  
  
    LIBS.translateX(JANGGUT_4_MODEL_MATRIX, LIBS.degToRad(25));
    LIBS.translateY(JANGGUT_4_MODEL_MATRIX, LIBS.degToRad(-6));
    LIBS.translateZ(JANGGUT_4_MODEL_MATRIX, LIBS.degToRad(53));
    LIBS.rotateX(JANGGUT_4_MODEL_MATRIX, LIBS.degToRad(-20));
  
  
    LIBS.translateX(JANGGUT_5_MODEL_MATRIX, LIBS.degToRad(-26));
    LIBS.translateY(JANGGUT_5_MODEL_MATRIX, LIBS.degToRad(-1));
    LIBS.translateZ(JANGGUT_5_MODEL_MATRIX, LIBS.degToRad(52));
    LIBS.rotateX(JANGGUT_5_MODEL_MATRIX, LIBS.degToRad(-20));
  
  
    LIBS.translateX(JANGGUT_6_MODEL_MATRIX, LIBS.degToRad(32));
    LIBS.translateY(JANGGUT_6_MODEL_MATRIX, LIBS.degToRad(0));
    LIBS.translateZ(JANGGUT_6_MODEL_MATRIX, LIBS.degToRad(49));
    LIBS.rotateX(JANGGUT_6_MODEL_MATRIX, LIBS.degToRad(-30));
  
  
    LIBS.translateX(JANGGUT_7_MODEL_MATRIX, LIBS.degToRad(-32));
    LIBS.translateY(JANGGUT_7_MODEL_MATRIX, LIBS.degToRad(-8));
    LIBS.translateZ(JANGGUT_7_MODEL_MATRIX, LIBS.degToRad(49));
    LIBS.rotateX(JANGGUT_7_MODEL_MATRIX, LIBS.degToRad(-20));
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING THE KALUNG_T.V ========================== */
    LIBS.translateX(TORUS_TV_MODEL_MATRIX, LIBS.degToRad(-63.5));
    LIBS.translateY(TORUS_TV_MODEL_MATRIX, LIBS.degToRad(1));
    LIBS.translateZ(TORUS_TV_MODEL_MATRIX, LIBS.degToRad(1));
    LIBS.rotateX(TORUS_TV_MODEL_MATRIX, LIBS.degToRad(-60));
    /*====================================================================== */
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  
  
  
  
  
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SETTING GUNTER $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    /*========================= SETTING BODY GUNTER ========================== */
    LIBS.rotateX(GUNTER_MODEL_MATRIX, 1.57);
    LIBS.translateX(GUNTER_MODEL_MATRIX, 2);
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING BODY GUNTER 2 ========================== */
    LIBS.translateY(GUNTER_MODEL_MATRIX2, 0.18); // Rotate around the Y-axis
    LIBS.translateZ(GUNTER_MODEL_MATRIX2, 0.2);
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING MOUTH GUNTER  ========================== */
    LIBS.translateZ(GUNTER_MOUTH_MODEL_MATRIX, -0.4); // Move the cylinder along the X-axis
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING RIGHT GUNTER EYE ========================== */
    LIBS.translateZ(GUNTER_EYE_MODEL_MATRIX, -0.6); // Rotate around the Y-axis
    LIBS.translateY(GUNTER_EYE_MODEL_MATRIX, 0.85);
    LIBS.translateX(GUNTER_EYE_MODEL_MATRIX, 0.3);
    LIBS.rotateX(GUNTER_EYE_MODEL_MATRIX, 3.2);
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING LEFT GUNTER EYE ========================== */
    LIBS.translateZ(GUNTER_EYE_MODEL_MATRIX2, -0.6); // Rotate around the Y-axis
    LIBS.translateY(GUNTER_EYE_MODEL_MATRIX2, 0.85);
    LIBS.translateX(GUNTER_EYE_MODEL_MATRIX2, -0.3);
    LIBS.rotateX(GUNTER_EYE_MODEL_MATRIX2, 3.2);
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING RIGHT GUNTER HAND ========================== */
    LIBS.translateX(GUNTER_HAND_MODEL_MATRIX, 0.8);
    LIBS.translateY(GUNTER_HAND_MODEL_MATRIX, 0.4);
    LIBS.rotateX(GUNTER_HAND_MODEL_MATRIX, -0.4);
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING LEFT GUNTER HAND ========================== */
    LIBS.translateX(GUNTER_HAND_MODEL_MATRIX2, -0.8);
    LIBS.translateY(GUNTER_HAND_MODEL_MATRIX2, 0.4);
    LIBS.rotateX(GUNTER_HAND_MODEL_MATRIX2, -0.4);
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING RIGHT GUNTER FEET ========================== */
    LIBS.translateZ(GUNTER_FEET_MODEL_MATRIX, 0.95);
    LIBS.translateY(GUNTER_FEET_MODEL_MATRIX, 0.9);
    LIBS.translateX(GUNTER_FEET_MODEL_MATRIX, 0.4);
    LIBS.rotateX(GUNTER_FEET_MODEL_MATRIX, -3);
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING LEFT GUNTER FEET ========================== */
    LIBS.translateZ(GUNTER_FEET_MODEL_MATRIX2, 0.95);
    LIBS.translateY(GUNTER_FEET_MODEL_MATRIX2, 0.9);
    LIBS.translateX(GUNTER_FEET_MODEL_MATRIX2, -0.4);
    LIBS.rotateX(GUNTER_FEET_MODEL_MATRIX2, -3);
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING TORUS GUNTER ========================== */
    LIBS.translateX(GUNTER_TORUS_MODEL_MATRIX, -1.4);
    LIBS.translateY(GUNTER_TORUS_MODEL_MATRIX, 0.05);
    /*========================================================================== */
    /*========================= SETTING PERMATA GUNTER ========================== */
    LIBS.translateY(GUNTER_PERMATA_MODEL_MATRIX, 1); // Move the cylinder along the X-axis
    LIBS.translateZ(GUNTER_PERMATA_MODEL_MATRIX, 0.2); // Move the cylinder along the X-axis
  
  
  
  
  
  
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SETTING JAKE $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    //=================== SETTING BODY CYLINDER ===========================/
    LIBS.translateX(MODEL_BODY_CYLINDER_MATRIX_JAKE, LIBS.degToRad(0)); // Move the cylinder along the X-axis
    LIBS.translateY(MODEL_BODY_CYLINDER_MATRIX_JAKE, LIBS.degToRad(-90)); // Move the cylinder along the Y-axis
    LIBS.translateZ(MODEL_BODY_CYLINDER_MATRIX_JAKE, LIBS.degToRad(0)); // Move the cylinder along the Z-axis
    LIBS.rotateX(MODEL_BODY_CYLINDER_MATRIX_JAKE, 0);
    //=====================================================================/
  
  
  
  
    //=================== SETTING BODY TOP ===========================/
    LIBS.translateX(MODEL_BODY_TOP_MATRIX_JAKE, LIBS.degToRad(0)); // Move the cylinder along the X-axis
    LIBS.translateY(MODEL_BODY_TOP_MATRIX_JAKE, LIBS.degToRad(-22)); // Move the cylinder along the Y-axis
    LIBS.translateZ(MODEL_BODY_TOP_MATRIX_JAKE, LIBS.degToRad(0)); // Move the cylinder along the Z-axis
    LIBS.rotateX(MODEL_BODY_TOP_MATRIX_JAKE, LIBS.degToRad(90));
    //=====================================================================/
  
  
  
  
    //=================== SETTING BODY BOTTOM ===========================/
    LIBS.translateX(MODEL_BODY_BOTTOM_MATRIX_JAKE, LIBS.degToRad(0)); // Move the cylinder along the X-axis
    LIBS.translateY(MODEL_BODY_BOTTOM_MATRIX_JAKE, LIBS.degToRad(-43)); // Move the cylinder along the Y-axis
    LIBS.translateZ(MODEL_BODY_BOTTOM_MATRIX_JAKE, LIBS.degToRad(0)); // Move the cylinder along the Z-axis
    LIBS.rotateX(MODEL_BODY_BOTTOM_MATRIX_JAKE, LIBS.degToRad(-90));
    //=====================================================================/
  
  
  
  
    //=================== SETTING LEFT ARM CYLINDER ===========================/
    LIBS.translateX(MODEL_LEFT_ARM_MATRIX_JAKE, LIBS.degToRad(-50)); // Move the cylinder along the X-axis
    LIBS.translateY(MODEL_LEFT_ARM_MATRIX_JAKE, LIBS.degToRad(-20)); // Move the cylinder along the Y-axis
    LIBS.translateZ(MODEL_LEFT_ARM_MATRIX_JAKE, LIBS.degToRad(0)); // Move the cylinder along the Z-axis
    LIBS.rotateZ(MODEL_LEFT_ARM_MATRIX_JAKE, 0.25);
    LIBS.rotateX(MODEL_LEFT_ARM_MATRIX_JAKE, 3.15);
    // LIBS.rotateY(MODEL_BODY_CYLINDER_MATRIX_JAKE, LIBS.degToRad((1 * dt) / 10))
    //=====================================================================/
  
  
    //=================== SETTING RIGHT ARM CYLINDER ===========================/
    LIBS.translateX(MODEL_RIGHT_ARM_MATRIX_JAKE, LIBS.degToRad(50)); // Move the cylinder along the X-axis
    LIBS.translateY(MODEL_RIGHT_ARM_MATRIX_JAKE, LIBS.degToRad(-20)); // Move the cylinder along the Y-axis
    LIBS.translateZ(MODEL_RIGHT_ARM_MATRIX_JAKE, LIBS.degToRad(0)); // Move the cylinder along the Z-axis
    LIBS.rotateZ(MODEL_RIGHT_ARM_MATRIX_JAKE, -0.25);
    LIBS.rotateX(MODEL_RIGHT_ARM_MATRIX_JAKE, 3.15);
    // LIBS.rotateY(MODEL_BODY_CYLINDER_MATRIX_JAKE, LIBS.degToRad((1 * dt) / 10))
    //=====================================================================/
  
  
    //=================== SETTING LEFT FEET CYLINDER ===========================/
    LIBS.translateX(MODEL_LEFT_FEET_MATRIX_JAKE, LIBS.degToRad(-20)); // Move the cylinder along the X-axis
    LIBS.translateY(MODEL_LEFT_FEET_MATRIX_JAKE, LIBS.degToRad(-50)); // Move the cylinder along the Y-axis
    LIBS.translateZ(MODEL_LEFT_FEET_MATRIX_JAKE, LIBS.degToRad(0)); // Move the cylinder along the Z-axis
    LIBS.rotateX(MODEL_LEFT_FEET_MATRIX_JAKE, LIBS.degToRad(180));
    // LIBS.rotateY(MODEL_BODY_CYLINDER_MATRIX_JAKE, LIBS.degToRad((1 * dt) / 10))
    //=====================================================================/
  
  
    //=================== SETTING RIGHT FEET CYLINDER ===========================/
    LIBS.translateX(MODEL_RIGHT_FEET_MATRIX_JAKE, LIBS.degToRad(20)); // Move the cylinder along the X-axis
    LIBS.translateY(MODEL_RIGHT_FEET_MATRIX_JAKE, LIBS.degToRad(-50)); // Move the cylinder along the Y-axis
    LIBS.translateZ(MODEL_RIGHT_FEET_MATRIX_JAKE, LIBS.degToRad(0)); // Move the cylinder along the Z-axis
    LIBS.rotateX(MODEL_RIGHT_FEET_MATRIX_JAKE, LIBS.degToRad(180));
    // LIBS.rotateY(MODEL_BODY_CYLINDER_MATRIX_JAKE, LIBS.degToRad((1 * dt) / 10))
    //=====================================================================/
  
  
  
  
    //=================== SETTING LEFT FOOT ===========================/
    LIBS.translateX(MODEL_LEFT_FOOT_MATRIX_JAKE, LIBS.degToRad(0)); // Move the cylinder along the X-axis
    LIBS.translateY(MODEL_LEFT_FOOT_MATRIX_JAKE, LIBS.degToRad(140)); // Move the cylinder along the Y-axis
    LIBS.translateZ(MODEL_LEFT_FOOT_MATRIX_JAKE, LIBS.degToRad(67)); // Move the cylinder along the Z-axis
    LIBS.rotateX(MODEL_LEFT_FOOT_MATRIX_JAKE, LIBS.degToRad(0));
    //=====================================================================/
  
  
  
  
    //=================== SETTING RIGHT FOOT ===========================/
    LIBS.translateX(MODEL_RIGHT_FOOT_MATRIX_JAKE, LIBS.degToRad(0)); // Move the cylinder along the X-axis
    LIBS.translateY(MODEL_RIGHT_FOOT_MATRIX_JAKE, LIBS.degToRad(140)); // Move the cylinder along the Y-axis
    LIBS.translateZ(MODEL_RIGHT_FOOT_MATRIX_JAKE, LIBS.degToRad(67)); // Move the cylinder along the Z-axis
    LIBS.rotateX(MODEL_RIGHT_FOOT_MATRIX_JAKE, LIBS.degToRad(0));
    //=====================================================================/
  
  
  
  
    // /*========================= SETTING THE LEFT EYE ========================== */
    LIBS.translateX(MODEL_LEFT_EYE_MATRIX_JAKE, LIBS.degToRad(-15)); // Move the cylinder along the X-axis
    LIBS.translateY(MODEL_LEFT_EYE_MATRIX_JAKE, LIBS.degToRad(0)); // Move the cylinder along the Y-axis
    LIBS.translateZ(MODEL_LEFT_EYE_MATRIX_JAKE, LIBS.degToRad(58)); // Move the cylinder along the Z-axis
    LIBS.rotateX(MODEL_LEFT_EYE_MATRIX_JAKE, LIBS.degToRad(90)); // Rotate the cylinder around the X-axis
    // /*====================================================================== */
  
  
  
  
    // /*========================= SETTING THE RIGHT EYE ========================== */
    LIBS.translateX(MODEL_RIGHT_EYE_MATRIX_JAKE, LIBS.degToRad(20)); // Move the cylinder along the X-axis
    LIBS.translateY(MODEL_RIGHT_EYE_MATRIX_JAKE, LIBS.degToRad(0)); // Move the cylinder along the Y-axis
    LIBS.translateZ(MODEL_RIGHT_EYE_MATRIX_JAKE, LIBS.degToRad(58)); // Move the cylinder along the Z-axis
    LIBS.rotateX(MODEL_RIGHT_EYE_MATRIX_JAKE, LIBS.degToRad(90)); // Rotate the cylinder around the X-axis
    // /*====================================================================== */
  
  
  
  
    /*========================= SETTING BLACK NOSE ========================== */
    LIBS.translateX(MODEL_BLACK_NOSE_MATRIX_JAKE, LIBS.degToRad(0));
    LIBS.translateY(MODEL_BLACK_NOSE_MATRIX_JAKE, LIBS.degToRad(-20));
    LIBS.translateZ(MODEL_BLACK_NOSE_MATRIX_JAKE, LIBS.degToRad(100));
    LIBS.rotateX(MODEL_BLACK_NOSE_MATRIX_JAKE, LIBS.degToRad(90));
    /*====================================================================== */
  
  
  
  
    /*========================= SETTING CURVE TORUS ========================== */
    LIBS.translateX(MODEL_CURVE_TORUS_MATRIX_JAKE, LIBS.degToRad(0));
    LIBS.translateY(MODEL_CURVE_TORUS_MATRIX_JAKE, LIBS.degToRad(-25));
    LIBS.translateZ(MODEL_CURVE_TORUS_MATRIX_JAKE, LIBS.degToRad(58));
    LIBS.rotateX(MODEL_CURVE_TORUS_MATRIX_JAKE, LIBS.degToRad(360));
    /*====================================================================== */
  
  
  
  
    //=================== SETTING LEFT EAR CYLINDER ===========================/
    LIBS.translateX(MODEL_LEFT_EAR_MATRIX_JAKE, LIBS.degToRad(0)); // Move the cylinder along the X-axis
    LIBS.translateY(MODEL_LEFT_EAR_MATRIX_JAKE, LIBS.degToRad(80)); // Move the cylinder along the Y-axis
    LIBS.translateZ(MODEL_LEFT_EAR_MATRIX_JAKE, LIBS.degToRad(0)); // Move the cylinder along the Z-axis
    LIBS.rotateY(MODEL_LEFT_EAR_MATRIX_JAKE, LIBS.degToRad(45));
    LIBS.rotateX(MODEL_LEFT_EAR_MATRIX_JAKE, LIBS.degToRad(-90));
    //=====================================================================/
  
  
  
  
    //=================== SETTING RIGHT EAR CYLINDER ===========================/
    LIBS.translateX(MODEL_RIGHT_EAR_MATRIX_JAKE, LIBS.degToRad(0)); // Move the cylinder along the X-axis
    LIBS.translateY(MODEL_RIGHT_EAR_MATRIX_JAKE, LIBS.degToRad(80)); // Move the cylinder along the Y-axis
    LIBS.translateZ(MODEL_RIGHT_EAR_MATRIX_JAKE, LIBS.degToRad(0)); // Move the cylinder along the Z-axis
    LIBS.rotateY(MODEL_RIGHT_EAR_MATRIX_JAKE, LIBS.degToRad(-45));
    LIBS.rotateX(MODEL_RIGHT_EAR_MATRIX_JAKE, LIBS.degToRad(-90));
    //=====================================================================/
  
  
  
  
    //============================= SETTING TAIL ===========================/
    LIBS.translateX(MODEL_TAIL_MATRIX_JAKE, LIBS.degToRad(0)); // Move the cylinder along the X-axis
    LIBS.translateY(MODEL_TAIL_MATRIX_JAKE, LIBS.degToRad(-70)); // Move the cylinder along the Y-axis
    LIBS.translateZ(MODEL_TAIL_MATRIX_JAKE, LIBS.degToRad(10)); // Move the cylinder along the Z-axis
    LIBS.rotateX(MODEL_TAIL_MATRIX_JAKE, LIBS.degToRad(-45))
    //=====================================================================/
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  
  
  
  
  
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ SETTING ENVIRONTMENT $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    //============================= SETTING BASE ===========================/
    LIBS.translateY(BASE_MODEL_MATRIX, -1);
    LIBS.translateZ(BASE_MODEL_MATRIX, 10);
    LIBS.rotateX(BASE_MODEL_MATRIX, 1.57);
    //=====================================================================/
  
  
    //============================= SETTING BASE 2 ===========================/
    LIBS.translateX(BASE2_MODEL_MATRIX, -15);
    LIBS.translateY(BASE2_MODEL_MATRIX, -1);
    LIBS.translateZ(BASE2_MODEL_MATRIX, -40);
    LIBS.rotateX(BASE2_MODEL_MATRIX, 1.57);
    //=====================================================================/
  
  
    //============================= SETTING BASE 3 ===========================/
    LIBS.translateX(BASE3_MODEL_MATRIX, LIBS.degToRad(-2000));
    LIBS.translateZ(BASE3_MODEL_MATRIX, LIBS.degToRad(-600));
    LIBS.rotateX(BASE3_MODEL_MATRIX, 1.57);
    //=====================================================================/
  
  
  
  
    //============================= SETTING BATANG BESAR ===========================/
    LIBS.translateY(BATANG_MODEL_MATRIX, -1);
    LIBS.translateX(BATANG_MODEL_MATRIX, 20);
    LIBS.rotateX(BATANG_MODEL_MATRIX, 11);
  
  
    LIBS.translateY(BATANG1_MODEL_MATRIX, -1);
    LIBS.translateX(BATANG1_MODEL_MATRIX, 25);
    LIBS.translateZ(BATANG1_MODEL_MATRIX, -1);
    LIBS.rotateZ(BATANG1_MODEL_MATRIX, 11);
    LIBS.rotateY(BATANG1_MODEL_MATRIX, 1);
  
  
    LIBS.translateY(BATANG2_MODEL_MATRIX, -1);
    LIBS.translateX(BATANG2_MODEL_MATRIX, 15);
    LIBS.translateZ(BATANG2_MODEL_MATRIX, -1);
    LIBS.rotateZ(BATANG2_MODEL_MATRIX, 11.5);
    LIBS.rotateY(BATANG2_MODEL_MATRIX, 2.3);
  
  
  
  
    LIBS.translateZ(BATANG3_MODEL_MATRIX, -6);
    LIBS.translateY(BATANG3_MODEL_MATRIX, -1);
    LIBS.translateX(BATANG3_MODEL_MATRIX, 30);
    LIBS.rotateX(BATANG3_MODEL_MATRIX, 11);
    LIBS.rotateY(BATANG3_MODEL_MATRIX, -0.8);
    //=====================================================================/
  
  
  
  
    //============================= SETTING DAUN POHON BESAR ===========================/
    LIBS.translateZ(DAUN_MODEL_MATRIX, -24);
    LIBS.translateX(DAUN_MODEL_MATRIX, 20);
  
  
    LIBS.translateZ(DAUN1_MODEL_MATRIX, -10);
    LIBS.translateY(DAUN1_MODEL_MATRIX, -1);
    LIBS.translateX(DAUN1_MODEL_MATRIX, 28)
  
  
    LIBS.translateZ(DAUN2_MODEL_MATRIX, -12);
    LIBS.translateY(DAUN2_MODEL_MATRIX, 3);
    LIBS.translateX(DAUN2_MODEL_MATRIX, 10);
  
  
    LIBS.translateZ(DAUN3_MODEL_MATRIX, -16);
    LIBS.translateY(DAUN3_MODEL_MATRIX, 0);
    LIBS.translateX(DAUN3_MODEL_MATRIX, 32);
    //=====================================================================/
    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  
  
  
    // LOKASI TUBUH JAKE:
    LIBS.translateX(MODEL_MATRIX_JAKE, LIBS.degToRad(-1400));
    LIBS.translateY(MODEL_MATRIX_JAKE, LIBS.degToRad(155));
    LIBS.translateZ(MODEL_MATRIX_JAKE, LIBS.degToRad(1040));
    LIBS.rotateY(MODEL_MATRIX_JAKE, LIBS.degToRad(-180));
  
  
    // LOKASI TUBUH GUNTER:
    LIBS.translateX(GUNTER_MODEL_MATRIX, LIBS.degToRad(290));
    LIBS.translateY(GUNTER_MODEL_MATRIX, LIBS.degToRad(-10));
    LIBS.translateZ(GUNTER_MODEL_MATRIX, LIBS.degToRad(-60));
    LIBS.rotateY(GUNTER_MODEL_MATRIX, LIBS.degToRad(-90));
  
  
    // LOKASI TUBUH TV:
    LIBS.translateX(MODEL_MATRIX, LIBS.degToRad(-900));
    LIBS.translateX(MODEL_MATRIX, LIBS.degToRad(-360));
    LIBS.translateY(MODEL_MATRIX, LIBS.degToRad(20));
    LIBS.rotateY(MODEL_MATRIX, LIBS.degToRad(90));
  
  
    // LOKASI Kamera
  
  
    LIBS.translateX(VIEW_MATRIX, 20);
    LIBS.translateY(VIEW_MATRIX, -3);
    LIBS.translateZ(VIEW_MATRIX, -20);
    
    // LIBS.translateY(VIEW_MATRIX, -4);
    // LIBS.translateZ(VIEW_MATRIX, -20);
    // LIBS.translateX(VIEW_MATRIX, 0);
    // LIBS.rotateX(VIEW_MATRIX, 0);
  
  
    // LOKASI LEMARI:
    LIBS.translateX(lemari_MODEL_MATRIX, LIBS.degToRad(-200));
    LIBS.translateY(lemari_MODEL_MATRIX, LIBS.degToRad(180));
    LIBS.translateZ(lemari_MODEL_MATRIX, 0);
    LIBS.rotateY(lemari_MODEL_MATRIX, LIBS.degToRad(-90));
  
  
    // LOKASI TANGGA:
    LIBS.rotateY(kiri_tangga_MODEL_MATRIX, LIBS.degToRad(88));
    LIBS.translateX(kiri_tangga_MODEL_MATRIX, LIBS.degToRad(-108));
    LIBS.translateZ(kiri_tangga_MODEL_MATRIX, 0);
    LIBS.translateY(kiri_tangga_MODEL_MATRIX, LIBS.degToRad(200));
    LIBS.rotateX(kiri_tangga_MODEL_MATRIX, LIBS.degToRad(0));
    LIBS.rotateZ(kiri_tangga_MODEL_MATRIX, LIBS.degToRad(10));
  
  
    //lokasi path
    LIBS.translateX(ZIG_MODEL_MATRIX, -10);
    LIBS.translateY(ZIG_MODEL_MATRIX, -20);
    LIBS.translateZ(ZIG_MODEL_MATRIX, 0);
    LIBS.rotateX(ZIG_MODEL_MATRIX, LIBS.degToRad(200));
    LIBS.rotateZ(ZIG_MODEL_MATRIX, LIBS.degToRad(-15));
  
  
    //lokasi mount
    LIBS.translateY(MOUNT_MODEL_MATRIX, -50);
    LIBS.translateX(MOUNT_MODEL_MATRIX, -19);
    LIBS.rotateX(MOUNT_MODEL_MATRIX, LIBS.degToRad(-90));
  
  
    LIBS.translateY(MOUNT_MODEL_MATRIX1, -52);
    LIBS.translateX(MOUNT_MODEL_MATRIX1, -5);
    LIBS.rotateX(MOUNT_MODEL_MATRIX1, LIBS.degToRad(-90));
  
  
    //lokaasi gua
    LIBS.translateZ(GUA_MODEL_MATRIX, -18);
    LIBS.translateY(GUA_MODEL_MATRIX, -40.5);
    LIBS.translateX(GUA_MODEL_MATRIX, -16.5);
    LIBS.rotateX(GUA_MODEL_MATRIX, LIBS.degToRad(-90));
  
  
    //lokasi mata ice kingdom
    LIBS.translateZ(EYE_KINGDOM_MODEL_MATRIX, -24);
    LIBS.translateY(EYE_KINGDOM_MODEL_MATRIX, -40.5);
    LIBS.translateX(EYE_KINGDOM_MODEL_MATRIX, -19.5);
    
    LIBS.translateZ(EYE_KINGDOM_MODEL_MATRIX1, -24);
    LIBS.translateY(EYE_KINGDOM_MODEL_MATRIX1, -42);
    LIBS.translateX(EYE_KINGDOM_MODEL_MATRIX1, -14);
  
  
    
  
  
    
  
  
  
  
    // LOKASI KALENG:
  
  
    GL.clearColor(107 / 255, 207 / 255, 240 / 255, 1);
    GL.enable(GL.DEPTH_TEST);
    GL.depthFunc(GL.LEQUAL);
  
  
    var timeRef = 0;
  
  
    var animate = function (time) {
        GL.viewport(0, 0, CANVAS.width, CANVAS.height);
        GL.clear(GL.COLOR_BUFFER_BIT);
  
  
        var dt = time - timeRef;
        timeRef = time;
  
  
  
  
        //######################################################################################################
        //############################ DRAWING T.V #############################################################
        //######################################################################################################
  
  
        /*========================= DRAWING THE SPHERE ========================= */
        GL.bindBuffer(GL.ARRAY_BUFFER, TRIANGLE_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 4 * (3 + 3), 0);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 4 * (3 + 3), 3 * 4);
        GL.uniformMatrix4fv(_PMatrix, false, PROJECTION_MATRIX);
        GL.uniformMatrix4fv(_VMatrix, false, VIEW_MATRIX);
        GL.uniformMatrix4fv(_MMatrix, false, MODEL_MATRIX);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, TRIANGLE_FACES);
        GL.enableVertexAttribArray(_color);
        GL.bindBuffer(GL.ARRAY_BUFFER, TRIANGLE_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLES, sphere.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING THE CONE ========================== */
        GL.bindBuffer(GL.ARRAY_BUFFER, CONE_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, CONE_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLES, 0, cone.vertices.length / 3);
        GL.flush();
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING THE LEFT EYE ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(CYLINDER_MODEL_MATRIX, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, EYE_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, EYE_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, eye.vertices.length / 3);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING THE RIGHT EYE ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(CYLINDER_MODEL_MATRIX_2, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, SECOND_EYE_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, SECOND_EYE_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, second_eye.vertices.length / 3);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING THE RIGHT BROW ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(Right_BROW_MODEL_MATRIX, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_BROW_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_BROW_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, right_brow.vertices.length / 3);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING THE LEFT BROW ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(Left_BROW_MODEL_MATRIX, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_BROW_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_BROW_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, left_brow.vertices.length / 3);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING THE LEFT HAND ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(Left_HAND_MODEL_MATRIX, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_HAND_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, LEFT_HAND_INDICES); // Bind indices buffer
        GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_HAND_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLE_STRIP, left_hand.indices.length, GL.UNSIGNED_SHORT, 0); // Use drawElements with GL.TRIANGLE_STRIP
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING THE RIGHT HAND ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(Right_HAND_MODEL_MATRIX, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_HAND_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, RIGHT_HAND_INDICES); // Bind indices buffer
        GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_HAND_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLE_STRIP, right_hand.indices.length, GL.UNSIGNED_SHORT, 0); // Use drawElements with GL.TRIANGLE_STRIP
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING THE LEFT FOOT ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(Left_FOOT_MODEL_MATRIX, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_FOOT_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, LEFT_FOOT_INDICES); // Bind indices buffer
        GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_FOOT_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLE_STRIP, left_foot.indices.length, GL.UNSIGNED_SHORT, 0); // Use drawElements with GL.TRIANGLE_STRIP
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING THE RIGHT FOOT ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(Right_FOOT_MODEL_MATRIX, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_FOOT_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, RIGHT_FOOT_INDICES); // Bind indices buffer
        GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_FOOT_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLE_STRIP, right_foot.indices.length, GL.UNSIGNED_SHORT, 0); // Use drawElements with GL.TRIANGLE_STRIP
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING THE TAIL ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(TAIL_MODEL_MATRIX, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, TAIL_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, TAIL_INDICES); // Bind indices buffer
        GL.bindBuffer(GL.ARRAY_BUFFER, TAIL_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLE_STRIP, tail.indices.length, GL.UNSIGNED_SHORT, 0); // Use drawElements with GL.TRIANGLE_STRIP
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING THE NOSE ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(NOSE_MODEL_MATRIX, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, NOSE_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, NOSE_INDICES); // Bind indices buffer
        GL.bindBuffer(GL.ARRAY_BUFFER, NOSE_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLE_STRIP, nose.indices.length, GL.UNSIGNED_SHORT, 0); // Use drawElements with GL.TRIANGLE_STRIP
        /*======================================================================= */
  
  
  
  
        /*========================= DRAWING THE NOSE 2 ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(NOSE_MODEL_MATRIX_2, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, NOSE_VERTEX_2);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, NOSE_INDICES_2); // Bind indices buffer
        GL.bindBuffer(GL.ARRAY_BUFFER, NOSE_COLOR_2);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLE_STRIP, nose_2.indices.length, GL.UNSIGNED_SHORT, 0); // Use drawElements with GL.TRIANGLE_STRIP
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING THE LEFT EAR ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(LEFT_EAR_MODEL_MATRIX, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_EAR_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, LEFT_EAR_INDICES); // Bind indices buffer
        GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_EAR_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLE_STRIP, left_ear.indices.length, GL.UNSIGNED_SHORT, 0); // Use drawElements with GL.TRIANGLE_STRIP
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING THE RIGHT EAR ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(RIGHT_EAR_MODEL_MATRIX, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_EAR_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, RIGHT_EAR_INDICES); // Bind indices buffer
        GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_EAR_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLE_STRIP, right_ear.indices.length, GL.UNSIGNED_SHORT, 0); // Use drawElements with GL.TRIANGLE_STRIP
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING THE JANGGUT ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(JANGGUT_1_MODEL_MATRIX, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_1_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_1_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, janggut_1.vertices.length / 3);
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(JANGGUT_2_MODEL_MATRIX, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_2_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_2_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, janggut_2.vertices.length / 3);
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(JANGGUT_3_MODEL_MATRIX, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_3_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_3_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, janggut_3.vertices.length / 3);
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(JANGGUT_4_MODEL_MATRIX, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_4_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_4_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, janggut_4.vertices.length / 3);
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(JANGGUT_5_MODEL_MATRIX, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_5_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_5_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, janggut_5.vertices.length / 3);
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(JANGGUT_6_MODEL_MATRIX, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_6_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_6_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, janggut_6.vertices.length / 3);
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(JANGGUT_7_MODEL_MATRIX, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_7_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, JANGGUT_7_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, janggut_7.vertices.length / 3);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING THE KALUNG ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(TORUS_TV_MODEL_MATRIX, MODEL_MATRIX));
        GL.bindBuffer(GL.ARRAY_BUFFER, TORUS_TV_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, TORUS_TV_INDICES);
        GL.bindBuffer(GL.ARRAY_BUFFER, TORUS_TV_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLES, torus_TV.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING THE MOUTH ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MULUT_MODEL_MATRIX, MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, MULUT_TV_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, MULUT_TV_INDICES); // Bind indices buffer
        GL.bindBuffer(GL.ARRAY_BUFFER, MULUT_TV_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLE_STRIP, torus_mulut_tv.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
        //######################################################################################################
        //######################################################################################################
        //######################################################################################################
  
  
  
  
  
  
        //######################################################################################################
        //############################ DRAWING GUNTER ##########################################################
        //######################################################################################################
  
  
        /*========================= DRAWING BODY GUNTER ========================== */
        GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_BODY_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 4 * (3 + 3), 0);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 4 * (3 + 3), 3 * 4);
        GL.uniformMatrix4fv(_PMatrix, false, PROJECTION_MATRIX);
        GL.uniformMatrix4fv(_VMatrix, false, VIEW_MATRIX);
        GL.uniformMatrix4fv(_MMatrix, false, GUNTER_MODEL_MATRIX); // Pass the new model matrix
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, GUNTER_BODY_FACES);
        GL.enableVertexAttribArray(_color);
        GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_BODY_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLES, gunter_body.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING BODY GUNTER 2 ========================== */
        GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_BODY_VERTEX2);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 4 * (3 + 3), 0);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 4 * (3 + 3), 3 * 4);
        GL.uniformMatrix4fv(_PMatrix, false, PROJECTION_MATRIX);
        GL.uniformMatrix4fv(_VMatrix, false, VIEW_MATRIX);
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(GUNTER_MODEL_MATRIX2, GUNTER_MODEL_MATRIX)); // Pass the new model matrix
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, GUNTER_BODY_FACES2);
        GL.enableVertexAttribArray(_color);
        GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_BODY_COLOR2);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLES, gunter_body2.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING MOUTH GUNTER ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(GUNTER_MOUTH_MODEL_MATRIX, GUNTER_MODEL_MATRIX));
        // GL.uniformMatrix4fv(_MMatrix, false, GUNTER_MOUTH_MODEL_MATRIX);
        GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_MOUTH_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_MOUTH_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLES, 0, GUNTER_MOUTH.vertices.length / 3);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING RIGHT GUNTER EYE ========================== */
        GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_EYE_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 4 * (3 + 3), 0);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 4 * (3 + 3), 3 * 4);
        GL.uniformMatrix4fv(_PMatrix, false, PROJECTION_MATRIX);
        GL.uniformMatrix4fv(_VMatrix, false, VIEW_MATRIX);
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(GUNTER_EYE_MODEL_MATRIX, GUNTER_MODEL_MATRIX)); // Pass the new model matrix
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, GUNTER_EYE_FACES);
        GL.enableVertexAttribArray(_color);
        GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_EYE_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLES, GUNTER_EYE.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING LEFT GUNTER EYE ========================== */
        GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_EYE_VERTEX2);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 4 * (3 + 3), 0);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 4 * (3 + 3), 3 * 4);
        GL.uniformMatrix4fv(_PMatrix, false, PROJECTION_MATRIX);
        GL.uniformMatrix4fv(_VMatrix, false, VIEW_MATRIX);
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(GUNTER_EYE_MODEL_MATRIX2, GUNTER_MODEL_MATRIX)); // Pass the new model matrix
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, GUNTER_EYE_FACES2);
        GL.enableVertexAttribArray(_color);
        GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_EYE_COLOR2);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLES, GUNTER_EYE.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING THE RIGHT GUNTER_HAND ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(GUNTER_HAND_MODEL_MATRIX, GUNTER_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_HAND_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, GUNTER_HAND_INDICES); // Bind indices buffer
        GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_HAND_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLE_STRIP, GUNTER_HAND.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING THE LEFT GUNTER_HAND ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(GUNTER_HAND_MODEL_MATRIX2, GUNTER_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_HAND_VERTEX2);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, GUNTER_HAND_INDICES2); // Bind indices buffer
        GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_HAND_COLOR2);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLE_STRIP, GUNTER_HAND.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING THE RIGHT GUNTER_FEET ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(GUNTER_FEET_MODEL_MATRIX, GUNTER_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_FEET_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, GUNTER_FEET_INDICES); // Bind indices buffer
        GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_FEET_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLE_STRIP, GUNTER_FEET.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING THE LEFT GUNTER_FEET ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(GUNTER_FEET_MODEL_MATRIX2, GUNTER_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_FEET_VERTEX2);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, GUNTER_FEET_INDICES2); // Bind indices buffer
        GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_FEET_COLOR2);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLE_STRIP, GUNTER_FEET.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING TORUS GUNTER ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(GUNTER_TORUS_MODEL_MATRIX, GUNTER_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_TORUS_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, GUNTER_TORUS_INDICES); // Bind indices buffer
        GL.bindBuffer(GL.ARRAY_BUFFER, GUNTER_TORUS_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLE_STRIP, GUNTER_TORUS.indices.length, GL.UNSIGNED_SHORT, 0);
        /*=================================================================== */
        /*========================= DRAWING PERMATA GUNTER ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(GUNTER_PERMATA_MODEL_MATRIX, GUNTER_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, PERMATA_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, PERMATA_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, permata.vertices.length / 3);
        //######################################################################################################
        //######################################################################################################
        //######################################################################################################
  
  
  
  
  
  
        //######################################################################################################
        //############################ DRAWING JAKE ############################################################
        //######################################################################################################
  
  
        /*========================= DRAWING BODY CYLINDER ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_BODY_CYLINDER_MATRIX_JAKE, MODEL_MATRIX_JAKE)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, BODY_CYLINDER_VERTEX_JAKE);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, BODY_CYLINDER_COLOR_JAKE);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, bodyCylinder_jake.vertices.length / 3);
        // LIBS.rotateY(MODEL_BODY_CYLINDER_MATRIX_JAKE, LIBS.degToRad((1 * dt) / 10));
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING BODY TOP ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_BODY_TOP_MATRIX_JAKE, MODEL_MATRIX_JAKE)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, BODY_TOP_VERTEX_JAKE);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, BODY_TOP_COLOR_JAKE);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, BODY_TOP_INDICES_JAKE);
        GL.drawElements(GL.TRIANGLES, bodytop_jake.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING BODY BOTTOM ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_BODY_BOTTOM_MATRIX_JAKE, MODEL_MATRIX_JAKE)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, BODY_BOTTOM_VERTEX_JAKE);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, BODY_BOTTOM_COLOR_JAKE);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, BODY_BOTTOM_INDICES_JAKE);
        GL.drawElements(GL.TRIANGLES, bodybottom_jake.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING LEFT ARM ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_LEFT_ARM_MATRIX_JAKE, MODEL_MATRIX_JAKE)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_ARM_CYLINDER_VERTEX_JAKE);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_ARM_CYLINDER_COLOR_JAKE);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, leftarmcylinder_jake.vertices.length / 3);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING RIGHT ARM ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_RIGHT_ARM_MATRIX_JAKE, MODEL_MATRIX_JAKE)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_ARM_CYLINDER_VERTEX_JAKE);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_ARM_CYLINDER_COLOR_JAKE);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, rightarmcylinder_jake.vertices.length / 3);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING LEFT FEET ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_LEFT_FEET_MATRIX_JAKE, MODEL_MATRIX_JAKE)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_FEET_CYLINDER_VERTEX_JAKE);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_FEET_CYLINDER_COLOR_JAKE);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, leftfeetcylinder_jake.vertices.length / 3);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING RIGHT FEET ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_RIGHT_FEET_MATRIX_JAKE, MODEL_MATRIX_JAKE)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_FEET_CYLINDER_VERTEX_JAKE);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_FEET_CYLINDER_COLOR_JAKE);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, rightfeetcylinder_jake.vertices.length / 3);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING LEFT FOOT ========================== */
        var temp = LIBS.multiply(MODEL_LEFT_FOOT_MATRIX_JAKE, MODEL_LEFT_FEET_MATRIX_JAKE)
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_LEFT_FOOT_MATRIX_JAKE, MODEL_LEFT_FEET_MATRIX_JAKE)); // Apply the cylinder's model matrix
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(temp, MODEL_MATRIX_JAKE));
        GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_FOOT_VERTEX_JAKE);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_FOOT_COLOR_JAKE);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, LEFT_FOOT_INDICES_JAKE);
        GL.drawElements(GL.TRIANGLES, leftfoot_jake.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING RIGHT FOOT ========================== */
        var temp = LIBS.multiply(MODEL_RIGHT_FOOT_MATRIX_JAKE, MODEL_RIGHT_FEET_MATRIX_JAKE)
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_RIGHT_FOOT_MATRIX_JAKE, MODEL_RIGHT_FEET_MATRIX_JAKE)); // Apply the cylinder's model matrix
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(temp, MODEL_MATRIX_JAKE));
        GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_FOOT_VERTEX_JAKE);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_FOOT_COLOR_JAKE);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, RIGHT_FOOT_INDICES_JAKE);
        GL.drawElements(GL.TRIANGLES, rightfoot_jake.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
        // /*========================= DRAWING THE LEFT EYE ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_LEFT_EYE_MATRIX_JAKE, MODEL_MATRIX_JAKE)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_EYE_VERTEX_JAKE);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_EYE_COLOR_JAKE);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, lefteye_jake.vertices.length / 3);
        // /*====================================================================== */
  
  
  
  
        // /*========================= DRAWING THE LEFT EYE ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_RIGHT_EYE_MATRIX_JAKE, MODEL_MATRIX_JAKE)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_EYE_VERTEX_JAKE);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_EYE_COLOR_JAKE);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, righteye_jake.vertices.length / 3);
        // /*====================================================================== */
  
  
  
  
        /*/*========================= DRAWING BLACK NOSE ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_BLACK_NOSE_MATRIX_JAKE, MODEL_MATRIX_JAKE)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, BLACK_NOSE_VERTEX_JAKE);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, BLACK_NOSE_INDICES_JAKE); // Bind indices buffer
        GL.bindBuffer(GL.ARRAY_BUFFER, BLACK_NOSE_COLOR_JAKE);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLE_STRIP, blacknose_jake.indices.length, GL.UNSIGNED_SHORT, 0); // Use drawElements with GL.TRIANGLE_STRIP
        // ======================================================================= //
  
  
  
  
        /*========================= DRAWING CURVE TORUS ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_CURVE_TORUS_MATRIX_JAKE, MODEL_MATRIX_JAKE)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, UNOSE_TORUS_VERTEX_JAKE);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, UNOSE_TORUS_INDICES_JAKE); // Bind indices buffer
        GL.bindBuffer(GL.ARRAY_BUFFER, UNOSE_TORUS_COLOR_JAKE);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLE_STRIP, unose_torus_jake.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING LEFT EAR ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_LEFT_EAR_MATRIX_JAKE, MODEL_MATRIX_JAKE)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_EAR_VERTEX_JAKE);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, LEFT_EAR_COLOR_JAKE);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, LEFT_EAR_INDICES_JAKE);
        GL.drawElements(GL.TRIANGLES, leftear_jake.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING RIGHT EAR ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_RIGHT_EAR_MATRIX_JAKE, MODEL_MATRIX_JAKE)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_EAR_VERTEX_JAKE);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, RIGHT_EAR_COLOR_JAKE);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, RIGHT_EAR_INDICES_JAKE);
        GL.drawElements(GL.TRIANGLES, rightear_jake.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING TAIL ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_TAIL_MATRIX_JAKE, MODEL_MATRIX_JAKE)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, TAIL_VERTEX_JAKE);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, TAIL_COLOR_JAKE);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, TAIL_INDICES_JAKE);
        GL.drawElements(GL.TRIANGLES, tail_jake.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
        //######################################################################################################
        //######################################################################################################
        //######################################################################################################
  
  
  
  
  
  
        //######################################################################################################
        //############################ DRAWING ENVIRONTMENT ####################################################
        //######################################################################################################
  
  
        /*========================= DRAWING THE BASE ========================= */
        GL.bindBuffer(GL.ARRAY_BUFFER, BASE_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 4 * (3 + 3), 0);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 4 * (3 + 3), 3 * 4);
        GL.uniformMatrix4fv(_PMatrix, false, PROJECTION_MATRIX);
        GL.uniformMatrix4fv(_VMatrix, false, VIEW_MATRIX);
        GL.uniformMatrix4fv(_MMatrix, false, BASE_MODEL_MATRIX);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, BASE_INDICES);
        GL.enableVertexAttribArray(_color);
        GL.bindBuffer(GL.ARRAY_BUFFER, BASE_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLES, base.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
        /*========================= DRAWING THE BASE 2 ========================= */
        GL.bindBuffer(GL.ARRAY_BUFFER, BASE2_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 4 * (3 + 3), 0);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 4 * (3 + 3), 3 * 4);
        GL.uniformMatrix4fv(_MMatrix, false, BASE2_MODEL_MATRIX);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, BASE2_INDICES);
        GL.enableVertexAttribArray(_color);
        GL.bindBuffer(GL.ARRAY_BUFFER, BASE2_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLES, base2.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
         /*========================= DRAWING THE BASE 3 ========================= */
         GL.bindBuffer(GL.ARRAY_BUFFER, BASE3_VERTEX);
         GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 4 * (3 + 3), 0);
         GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 4 * (3 + 3), 3 * 4);
         GL.uniformMatrix4fv(_MMatrix, false, BASE3_MODEL_MATRIX);
         GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, BASE3_INDICES);
         GL.enableVertexAttribArray(_color);
         GL.bindBuffer(GL.ARRAY_BUFFER, BASE3_COLOR);
         GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
         GL.drawElements(GL.TRIANGLES, base3.indices.length, GL.UNSIGNED_SHORT, 0);
         /*====================================================================== */
  
  
        /*========================= DRAWING BATANG ========================= */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(BATANG_MODEL_MATRIX, BASE_MODEL_MATRIX));
        // GL.uniformMatrix4fv(_MMatrix, false, GUNTER_MOUTH_MODEL_MATRIX);
        GL.bindBuffer(GL.ARRAY_BUFFER, BATANG_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, BATANG_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLES, 0, batang.vertices.length / 3);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING BATANG 1 ========================= */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(BATANG1_MODEL_MATRIX, BASE_MODEL_MATRIX));
        // GL.uniformMatrix4fv(_MMatrix, false, GUNTER_MOUTH_MODEL_MATRIX);
        GL.bindBuffer(GL.ARRAY_BUFFER, BATANG1_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, BATANG1_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLES, 0, batang1.vertices.length / 3);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING BATANG 2========================= */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(BATANG2_MODEL_MATRIX, BASE_MODEL_MATRIX));
        // GL.uniformMatrix4fv(_MMatrix, false, GUNTER_MOUTH_MODEL_MATRIX);
        GL.bindBuffer(GL.ARRAY_BUFFER, BATANG2_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, BATANG2_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLES, 0, batang2.vertices.length / 3);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING BATANG3 ========================= */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(BATANG3_MODEL_MATRIX, BASE_MODEL_MATRIX));
        // GL.uniformMatrix4fv(_MMatrix, false, GUNTER_MOUTH_MODEL_MATRIX);
        GL.bindBuffer(GL.ARRAY_BUFFER, BATANG3_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, BATANG3_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLES, 0, batang3.vertices.length / 3);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING DAUN ========================= */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(DAUN_MODEL_MATRIX, BASE_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, DAUN_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, DAUN_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, DAUN_INDICES);
        GL.drawElements(GL.TRIANGLES, daun.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING DAUN 1 ========================= */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(DAUN1_MODEL_MATRIX, BASE_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, DAUN1_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, DAUN1_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, DAUN1_INDICES);
        GL.drawElements(GL.TRIANGLES, daun1.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING DAUN 2 ========================= */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(DAUN2_MODEL_MATRIX, BASE_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, DAUN2_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, DAUN2_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, DAUN2_INDICES);
        GL.drawElements(GL.TRIANGLES, daun2.indices.length, GL.UNSIGNED_SHORT, 0);
        /*====================================================================== */
  
  
  
  
        /*========================= DRAWING DAUN 3 ========================= */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(DAUN3_MODEL_MATRIX, BASE_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, DAUN3_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, DAUN3_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, DAUN3_INDICES);
        GL.drawElements(GL.TRIANGLES, daun3.indices.length, GL.UNSIGNED_SHORT, 0);
        // ================================================================== //
  
  
  
  
        //=========================== POHON 1 ==================================
        //-------------------------- BATANG 1 ----------------------------------
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_BATANG_POHON1_MATRIX, BASE_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, BATANG_POHON1_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, BATANG_POHON1_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, batangpohon1.vertices.length / 3);
        //---------------------------- DAUN 1 ----------------------------------
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_DAUN_POHON1_MATRIX, BASE_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, DAUN_POHON1_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, DAUN_POHON1_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, DAUN_POHON1_FACES);
        GL.drawElements(GL.TRIANGLES, daunpohon1.indices.length, GL.UNSIGNED_SHORT, 0);
        //======================================================================
  //---------------------------- PATH ----------------------------------
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(ZIG_MODEL_MATRIX, BASE_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, ZIG_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, ZIG_INDICES); // Bind indices buffer
        GL.bindBuffer(GL.ARRAY_BUFFER, ZIG_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLE_STRIP, zig.indices.length, GL.UNSIGNED_SHORT, 0); 
  
  
                //======================================================================
  //---------------------------- MOUNT ----------------------------------
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MOUNT_MODEL_MATRIX, BASE_MODEL_MATRIX));
        // GL.uniformMatrix4fv(_MMatrix, false, GUNTER_MOUTH_MODEL_MATRIX);
        GL.bindBuffer(GL.ARRAY_BUFFER, MOUNT_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, MOUNT_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLES, 0, mount.vertices.length / 3);
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MOUNT_MODEL_MATRIX1, BASE_MODEL_MATRIX));
        // GL.uniformMatrix4fv(_MMatrix, false, GUNTER_MOUTH_MODEL_MATRIX);
        GL.bindBuffer(GL.ARRAY_BUFFER, MOUNT_VERTEX1);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, MOUNT_COLOR1);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLES, 0, mount1.vertices.length / 3);
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(GUA_MODEL_MATRIX, BASE_MODEL_MATRIX));
        // GL.uniformMatrix4fv(_MMatrix, false, GUNTER_MOUTH_MODEL_MATRIX);
        GL.bindBuffer(GL.ARRAY_BUFFER, GUA_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, GUA_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLES, 0, gua.vertices.length / 3);
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(EYE_KINGDOM_MODEL_MATRIX, BASE_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, EYE_KINGDOM_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, EYE_KINGDOM_INDICES); // Bind indices buffer
        GL.bindBuffer(GL.ARRAY_BUFFER, EYE_KINGDOM_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLE_STRIP, eye_kingdom.indices.length, GL.UNSIGNED_SHORT, 0); // Use drawElements with GL.TRIANGLE_STRIP
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(EYE_KINGDOM_MODEL_MATRIX1, BASE_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, EYE_KINGDOM_VERTEX1);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, EYE_KINGDOM_INDICES1); // Bind indices buffer
        GL.bindBuffer(GL.ARRAY_BUFFER, EYE_KINGDOM_COLOR1);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLE_STRIP, eye_kingdom.indices.length, GL.UNSIGNED_SHORT, 0); // Use drawElements with GL.TRIANGLE_STRIP
    
  
  
  
  
        //============================ KALENG 1 ======================================
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_KALENG1_MATRIX, BASE_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, KALENG1_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, KALENG1_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, kaleng1.vertices.length / 3);
        //============================================================================
  
  
        //============================ KALENG 2 ======================================
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_KALENG2_MATRIX, BASE_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, KALENG2_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, KALENG2_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, kaleng2.vertices.length / 3);
        //============================================================================
  
  
        //============================ MATAHARI ====================================
        GL.bindBuffer(GL.ARRAY_BUFFER, MATAHARI_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 4 * (3 + 3), 0);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 4 * (3 + 3), 3 * 4);
        GL.uniformMatrix4fv(_PMatrix, false, PROJECTION_MATRIX);
        GL.uniformMatrix4fv(_VMatrix, false, VIEW_MATRIX);
        GL.uniformMatrix4fv(_MMatrix, false, MODEL_MATAHARI_MATRIX);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, MATAHARI_FACES);
        GL.enableVertexAttribArray(_color);
        GL.bindBuffer(GL.ARRAY_BUFFER, MATAHARI_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawElements(GL.TRIANGLES, matahari.indices.length, GL.UNSIGNED_SHORT, 0);
        //============================================================================
  
  
  
  
  
  
        //================================= LEMARI & TANGGA ===================================================
        //----------------------------------- Draw lemari ------------------------------------------------
        GL.uniformMatrix4fv(_MMatrix, false, lemari_MODEL_MATRIX);
        GL.bindBuffer(GL.ARRAY_BUFFER, lemari_vertex);
        GL.vertexAttribPointer(position_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, lemari_color);
        GL.vertexAttribPointer(color_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, lemari_faces);
        GL.drawElements(GL.TRIANGLES, lemari.indices.length, GL.UNSIGNED_SHORT, 0);
        //---------------------------------=---------------------------------------------------------------
  
  
  
  
        //---------------------------------- Draw laci lemari ---------------------------------------------
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(laci_lemari_1_MODEL_MATRIX, lemari_MODEL_MATRIX));
        GL.bindBuffer(GL.ARRAY_BUFFER, laci_lemari_1_vertex);
        GL.vertexAttribPointer(position_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, laci_lemari_1_color);
        GL.vertexAttribPointer(color_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, laci_lemari_1_faces);
        GL.drawElements(GL.TRIANGLES, laci_lemari_1.indices.length, GL.UNSIGNED_SHORT, 0);
  
  
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(laci_lemari_2_MODEL_MATRIX, lemari_MODEL_MATRIX));
        GL.bindBuffer(GL.ARRAY_BUFFER, laci_lemari_2_vertex);
        GL.vertexAttribPointer(position_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, laci_lemari_2_color);
        GL.vertexAttribPointer(color_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, laci_lemari_2_faces);
        GL.drawElements(GL.TRIANGLES, laci_lemari_2.indices.length, GL.UNSIGNED_SHORT, 0);
        //----------------------------------------------------------------------------------------------------
  
  
  
  
        //------------------------------------- Draw kiri tangga ---------------------------------------------
        GL.uniformMatrix4fv(_MMatrix, false, kiri_tangga_MODEL_MATRIX);
        GL.bindBuffer(GL.ARRAY_BUFFER, kiri_tangga_vertex);
        GL.vertexAttribPointer(position_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, kiri_tangga_color);
        GL.vertexAttribPointer(color_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, kiri_tangga_faces);
        GL.drawElements(GL.TRIANGLES, kiri_tangga.indices.length, GL.UNSIGNED_SHORT, 0);
        //----------------------------------------------------------------------------------------------------
  
  
  
  
        //-------------------------------------- Draw tengah tangga ------------------------------------------
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(tengah_tangga_1_MODEL_MATRIX, kiri_tangga_MODEL_MATRIX));
        GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_1_vertex);
        GL.vertexAttribPointer(position_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_1_color);
        GL.vertexAttribPointer(color_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, tengah_tangga_1_faces);
        GL.drawElements(GL.TRIANGLES, tengah_tangga_1.indices.length, GL.UNSIGNED_SHORT, 0);
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(tengah_tangga_2_MODEL_MATRIX, kiri_tangga_MODEL_MATRIX));
        GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_2_vertex);
        GL.vertexAttribPointer(position_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_2_color);
        GL.vertexAttribPointer(color_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, tengah_tangga_2_faces);
        GL.drawElements(GL.TRIANGLES, tengah_tangga_2.indices.length, GL.UNSIGNED_SHORT, 0);
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(tengah_tangga_3_MODEL_MATRIX, kiri_tangga_MODEL_MATRIX));
        GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_3_vertex);
        GL.vertexAttribPointer(position_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_3_color);
        GL.vertexAttribPointer(color_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, tengah_tangga_3_faces);
        GL.drawElements(GL.TRIANGLES, tengah_tangga_3.indices.length, GL.UNSIGNED_SHORT, 0);
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(tengah_tangga_4_MODEL_MATRIX, kiri_tangga_MODEL_MATRIX));
        GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_4_vertex);
        GL.vertexAttribPointer(position_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_4_color);
        GL.vertexAttribPointer(color_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, tengah_tangga_4_faces);
        GL.drawElements(GL.TRIANGLES, tengah_tangga_4.indices.length, GL.UNSIGNED_SHORT, 0);
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(tengah_tangga_5_MODEL_MATRIX, kiri_tangga_MODEL_MATRIX));
        GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_5_vertex);
        GL.vertexAttribPointer(position_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_5_color);
        GL.vertexAttribPointer(color_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, tengah_tangga_5_faces);
        GL.drawElements(GL.TRIANGLES, tengah_tangga_5.indices.length, GL.UNSIGNED_SHORT, 0);
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(tengah_tangga_6_MODEL_MATRIX, kiri_tangga_MODEL_MATRIX));
        GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_6_vertex);
        GL.vertexAttribPointer(position_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_6_color);
        GL.vertexAttribPointer(color_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, tengah_tangga_6_faces);
        GL.drawElements(GL.TRIANGLES, tengah_tangga_6.indices.length, GL.UNSIGNED_SHORT, 0);
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(tengah_tangga_7_MODEL_MATRIX, kiri_tangga_MODEL_MATRIX));
        GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_7_vertex);
        GL.vertexAttribPointer(position_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_7_color);
        GL.vertexAttribPointer(color_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, tengah_tangga_7_faces);
        GL.drawElements(GL.TRIANGLES, tengah_tangga_7.indices.length, GL.UNSIGNED_SHORT, 0);
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(tengah_tangga_8_MODEL_MATRIX, kiri_tangga_MODEL_MATRIX));
        GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_8_vertex);
        GL.vertexAttribPointer(position_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_8_color);
        GL.vertexAttribPointer(color_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, tengah_tangga_8_faces);
        GL.drawElements(GL.TRIANGLES, tengah_tangga_8.indices.length, GL.UNSIGNED_SHORT, 0);
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(tengah_tangga_9_MODEL_MATRIX, kiri_tangga_MODEL_MATRIX));
        GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_9_vertex);
        GL.vertexAttribPointer(position_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_9_color);
        GL.vertexAttribPointer(color_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, tengah_tangga_9_faces);
        GL.drawElements(GL.TRIANGLES, tengah_tangga_9.indices.length, GL.UNSIGNED_SHORT, 0);
  
  
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(tengah_tangga_10_MODEL_MATRIX, kiri_tangga_MODEL_MATRIX));
        GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_10_vertex);
        GL.vertexAttribPointer(position_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, tengah_tangga_10_color);
        GL.vertexAttribPointer(color_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, tengah_tangga_10_faces);
        GL.drawElements(GL.TRIANGLES, tengah_tangga_10.indices.length, GL.UNSIGNED_SHORT, 0);
        //----------------------------------------------------------------------------------------------------
  
  
        //--------------------------------- Draw kanan tangga ------------------------------------------------
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(kanan_tangga_MODEL_MATRIX, kiri_tangga_MODEL_MATRIX));
        GL.bindBuffer(GL.ARRAY_BUFFER, kanan_tangga_vertex);
        GL.vertexAttribPointer(position_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, kanan_tangga_color);
        GL.vertexAttribPointer(color_vao, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, kanan_tangga_faces);
        GL.drawElements(GL.TRIANGLES, kanan_tangga.indices.length, GL.UNSIGNED_SHORT, 0);
        //------------------------------------------------------------------------------------------------------
        //============================================================================================
  
  
  
  
        //--------------------------------- KINCIR BAWAH --------------------------------------------------------
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_KINCIR_BAWAH_MATRIX, BASE_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, KINCIR_BAWAH_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, KINCIR_BAWAH_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, KINCIR_BAWAH_FACES);
        GL.drawElements(GL.TRIANGLES, kincirbawah.indices.length, GL.UNSIGNED_SHORT, 0);
        //-------------------------------------------------------------------------------------------------------
  
  
  
  
        //--------------------------------- KINCIR ATAS --------------------------------------------------------
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_KINCIR_ATAS_MATRIX, BASE_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, KINCIR_ATAS_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, KINCIR_ATAS_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ELEMENT_ARRAY_BUFFER, KINCIR_ATAS_FACES);
        GL.drawElements(GL.TRIANGLES, kinciratas.indices.length, GL.UNSIGNED_SHORT, 0);
        //-------------------------------------------------------------------------------------------------------
  
  
  
  
        //---------------------------------- TABUNG KINCIR -------------------------------------------------------
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_TABUNG_KINCIR_MATRIX, BASE_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, TABUNG_KINCIR_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, TABUNG_KINCIR_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, tabungkincir.vertices.length / 3);
        //-------------------------------------------------------------------------------------------------------
  
  
        /*========================= DRAWING KINCIR1 ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_KINCIR1_MATRIX, BASE_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, KINCIR1_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, KINCIR1_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, kincir1.vertices.length / 3);
        // /*====================================================================== */
  
  
        /*========================= DRAWING KINCIR2 ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_KINCIR2_MATRIX, BASE_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, KINCIR2_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, KINCIR2_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, kincir2.vertices.length / 3);
        /*====================================================================== */
  
  
        /*========================= DRAWING KINCIR3 ========================== */
        GL.uniformMatrix4fv(_MMatrix, false, LIBS.multiply(MODEL_KINCIR3_MATRIX, BASE_MODEL_MATRIX)); // Apply the cylinder's model matrix
        GL.bindBuffer(GL.ARRAY_BUFFER, KINCIR3_VERTEX);
        GL.vertexAttribPointer(_position, 3, GL.FLOAT, false, 0, 0);
        GL.bindBuffer(GL.ARRAY_BUFFER, KINCIR3_COLOR);
        GL.vertexAttribPointer(_color, 3, GL.FLOAT, false, 0, 0);
        GL.drawArrays(GL.TRIANGLE_STRIP, 0, kincir3.vertices.length / 3);
        /*====================================================================== */
  
  
  
  
        //######################################################################################################
        //######################################################################################################
        //###################################################################################################### 
  
  
        function TV_MOVE() {
          /*========================= ANIMATION Kaki berjalan========================= */
          if (time % 1000 < 500) {
              LIBS.rotateX(Left_FOOT_MODEL_MATRIX, 0.0090);
              LIBS.rotateX(Right_FOOT_MODEL_MATRIX, -0.0090);
          } else {
              LIBS.rotateX(Left_FOOT_MODEL_MATRIX, -0.0090);
              LIBS.rotateX(Right_FOOT_MODEL_MATRIX, 0.0090);
          }
          /*=========================================================================== */
  
  
  
  
          /*========================= ANIMATION tangan ========================= */
          if (time % 1000 < 500) {
              LIBS.rotateX(Left_HAND_MODEL_MATRIX, 0.0090);
              LIBS.rotateX(Right_HAND_MODEL_MATRIX, -0.0090);
          }
          else {
              LIBS.rotateX(Left_HAND_MODEL_MATRIX, -0.0090);
              LIBS.rotateX(Right_HAND_MODEL_MATRIX, 0.0090);
          }
          /*====================================================================== */
  
  
  
  
          /*========================= ANIMATION EARS ========================= */
          if (time % 600 < 300) {
              LIBS.rotateZ(LEFT_EAR_MODEL_MATRIX, 0.001);
              LIBS.rotateZ(RIGHT_EAR_MODEL_MATRIX, -0.001);
          }
          else {
              LIBS.rotateZ(LEFT_EAR_MODEL_MATRIX, -0.001);
              LIBS.rotateZ(RIGHT_EAR_MODEL_MATRIX, 0.001);
          }
          /*====================================================================== */
  
  
  
  
          /*========================= ANIMATION Tail ========================= */
          if (time % 1000 < 500) {
              LIBS.rotateX(TAIL_MODEL_MATRIX, 0.01);
          }
          else {
              LIBS.rotateX(TAIL_MODEL_MATRIX, -0.01);
          }
          /*====================================================================== */
        }
  
  
  
  
        function GUNTER_MOVE() {
          /*========================= ANIMATION EARS ========================= */
          if (time % 1000 < 500) {
            LIBS.rotateX(GUNTER_HAND_MODEL_MATRIX, 0.01);
            LIBS.rotateX(GUNTER_HAND_MODEL_MATRIX2, -0.01);
          }
          else {
            LIBS.rotateX(GUNTER_HAND_MODEL_MATRIX, -0.01);
            LIBS.rotateX(GUNTER_HAND_MODEL_MATRIX2, 0.01);
          }
  
  
          if (time % 1000 < 500) {
            LIBS.rotateZ(GUNTER_MODEL_MATRIX, 0.001);
            LIBS.rotateX(GUNTER_FEET_MODEL_MATRIX, -0.001);
            LIBS.rotateX(GUNTER_FEET_MODEL_MATRIX2, 0.001);
          }
          else {
            LIBS.rotateZ(GUNTER_MODEL_MATRIX, -0.001);
            LIBS.rotateX(GUNTER_FEET_MODEL_MATRIX, 0.001);
            LIBS.rotateX(GUNTER_FEET_MODEL_MATRIX2, -0.001);
          }
          /*====================================================================== */
        }
  
  
  
  
              function JAKE_MOVE() {
              /*========================= ANIMATION Kaki berjalan========================= */
              if (time % 1000 < 500) {
                  LIBS.rotateX(MODEL_LEFT_FEET_MATRIX_JAKE, -0.0050);
                  LIBS.rotateX(MODEL_RIGHT_FEET_MATRIX_JAKE, 0.0050);
  
  
  
  
              } else {
                  LIBS.rotateX(MODEL_LEFT_FEET_MATRIX_JAKE, 0.0050);
                  LIBS.rotateX(MODEL_RIGHT_FEET_MATRIX_JAKE, -0.0050);
  
  
  
  
              }
              /*=========================================================================== */
  
  
  
  
              /*========================= ANIMATION tangan ========================= */
              if (time % 1000 < 500) {
  
  
  
  
                  LIBS.rotateX(MODEL_LEFT_ARM_MATRIX_JAKE, 0.0050);
                  LIBS.rotateX(MODEL_RIGHT_ARM_MATRIX_JAKE, -0.0050);
              }
              else {
                  LIBS.rotateX(MODEL_LEFT_ARM_MATRIX_JAKE, -0.0050);
                  LIBS.rotateX(MODEL_RIGHT_ARM_MATRIX_JAKE, 0.0050);
              }
  
  
  
  
  
  
  
  
              /*========================= ANIMATION EARS ========================= */
              if (time % 600 < 300) {
                  LIBS.rotateZ(MODEL_LEFT_EAR_MATRIX_JAKE, 0.001);
                  LIBS.rotateZ(MODEL_RIGHT_EAR_MATRIX_JAKE, -0.001);
              }
              else {
                  LIBS.rotateZ(MODEL_LEFT_EAR_MATRIX_JAKE, -0.001);
                  LIBS.rotateZ(MODEL_RIGHT_EAR_MATRIX_JAKE, 0.001);
              }
  
  
  
  
              /*========================= ANIMATION Tail ========================= */
              if (time % 1000 < 500) {
                  LIBS.rotateX(MODEL_TAIL_MATRIX_JAKE, 0.01);
              }
              else {
                  LIBS.rotateX(MODEL_TAIL_MATRIX_JAKE, -0.01);
              }
  
  
  
  
          }
  
  
  
  
          // var radius = 50;
          // var pos_x = radius * Math.cos(time * 0.001);
          // var pos_y = radius * Math.sin(time * 0.001);
          // LIBS.translateX(MODEL_MATAHARI_MATRIX, pos_x, 20, pos_y);
  
  
  
  
          // TV
          if (time < 15100) {
              TV_MOVE();
              LIBS.translateX(MODEL_MATRIX, dt * 0.001);
              if (time < 7000) {
                  LIBS.translateX(VIEW_MATRIX, -dt * 0.001);
              }
          } else if (time > 15100 && time < 16800) {
              TV_MOVE();
              LIBS.rotateY(MODEL_MATRIX, -dt * 0.001);
          } else if (time > 16900 && time < 18200) {
              TV_MOVE();
              LIBS.translateZ(MODEL_MATRIX, dt * 0.001);
          } else if (time > 18200 && time < 19900) {
              TV_MOVE();
              LIBS.rotateY(MODEL_MATRIX, dt * 0.001);
          } else if (time > 19900 && time < 22100) {
              LIBS.rotateX(Left_HAND_MODEL_MATRIX, -dt * 0.001);
              LIBS.rotateX(Right_HAND_MODEL_MATRIX, -dt * 0.001);
              if (time > 20900 && time < 21400) {
                  LIBS.translateY(Left_HAND_MODEL_MATRIX, -dt * 0.001);
                  LIBS.translateY(Right_HAND_MODEL_MATRIX, -dt * 0.001);
              }
          } else if (time > 22100 && time < 22900) {
              LIBS.translateY(MODEL_MATRIX, dt * 0.001);
          } else if (time > 22900 && time < 23700) {
              LIBS.translateY(MODEL_MATRIX, -dt * 0.001);
          } else if (time > 23700 && time < 24500) {
              LIBS.translateY(MODEL_MATRIX, dt * 0.001);
          } else if (time > 24500 && time < 25300) {
              LIBS.translateY(MODEL_MATRIX, -dt * 0.001);
          } else if (time > 25300 && time < 26100) {
              LIBS.translateY(MODEL_MATRIX, dt * 0.001);
          } else if (time > 26100 && time < 26900) {
              LIBS.translateY(MODEL_MATRIX, -dt * 0.001);
          } else if (time > 26900 && time < 27700) {
              LIBS.translateY(MODEL_MATRIX, dt * 0.001);
          } else if (time > 27700 && time < 28500) {
              LIBS.translateY(MODEL_MATRIX, -dt * 0.001);
          } else if (time > 28500 && time < 29300) {
              LIBS.translateY(MODEL_MATRIX, dt * 0.001);
          } else if (time > 29300 && time < 30100) {
              LIBS.translateY(MODEL_MATRIX, -dt * 0.001);
          } else if (time > 30100 && time < 30900) {
              LIBS.translateY(MODEL_MATRIX, dt * 0.001);
          } else if (time > 30900 && time < 31700) {
              LIBS.translateY(MODEL_MATRIX, -dt * 0.001);
          } else if (time > 31700 && time < 32500) {
              LIBS.translateY(MODEL_MATRIX, dt * 0.001);
          } else if (time > 32500 && time < 33300) {
              LIBS.translateY(MODEL_MATRIX, -dt * 0.001);
          } else if (time > 33300 && time < 34100) {
              LIBS.translateY(MODEL_MATRIX, dt * 0.001);
          } else if (time > 34100 && time < 34900) {
              LIBS.translateY(MODEL_MATRIX, -dt * 0.001);
          } else if (time > 34900 && time < 35700) {
              LIBS.translateY(MODEL_MATRIX, dt * 0.001);
          } else if (time > 35700 && time < 36500) {
              LIBS.translateY(MODEL_MATRIX, -dt * 0.001);
          } else if (time > 36500 && time < 37300) {
              LIBS.translateY(MODEL_MATRIX, dt * 0.001);
          } else if (time > 37300 && time < 38100) {
              LIBS.translateY(MODEL_MATRIX, -dt * 0.001);
          } else if (time > 38100 && time < 38900) {
              LIBS.translateY(MODEL_MATRIX, dt * 0.001);
          } else if (time > 38900 && time < 39700) {
              LIBS.translateY(MODEL_MATRIX, -dt * 0.001);
          } else if (time > 39700 && time < 40500) {
              LIBS.translateY(MODEL_MATRIX, dt * 0.001);
          } else if (time > 40500 && time < 41300) {
              LIBS.translateY(MODEL_MATRIX, -dt * 0.001);
          } else if (time > 41300 && time < 42100) {
              LIBS.translateY(MODEL_MATRIX, dt * 0.001);
          } else if (time > 42100 && time < 42900) {
              LIBS.translateY(MODEL_MATRIX, -dt * 0.001);
          } else if (time > 42900 && time < 43700) {
              LIBS.translateY(MODEL_MATRIX, dt * 0.001);
          } else if (time > 43700 && time < 44500) {
              LIBS.translateY(MODEL_MATRIX, -dt * 0.001);
          } else if (time > 44500 && time < 46700) {
              LIBS.rotateX(Left_HAND_MODEL_MATRIX, dt * 0.001);
              LIBS.rotateX(Right_HAND_MODEL_MATRIX, dt * 0.001);
              if (time > 45500 && time < 46000) {
                  LIBS.translateY(Left_HAND_MODEL_MATRIX, dt * 0.001);
                  LIBS.translateY(Right_HAND_MODEL_MATRIX, dt * 0.001);
              }
          } else if (time > 47500 && time < 50800) {
              TV_MOVE();
              LIBS.rotateY(MODEL_MATRIX, -dt * 0.001);
          } else if (time > 51200 && time < 52600) {
              TV_MOVE();
              LIBS.rotateY(MODEL_MATRIX, -dt * 0.001);
          } else if (time > 52600 && time < 54300) {
              TV_MOVE();
              LIBS.translateZ(MODEL_MATRIX, -dt * 0.001);
          } else if (time > 54300 && time < 57800) {
              TV_MOVE();
              LIBS.rotateY(MODEL_MATRIX, -dt * 0.001);
          } else if (time > 57800 && time < 58300) {
              TV_MOVE();
              LIBS.rotateY(MODEL_MATRIX, dt * 0.0004);
          }
  
  
  
  
          if (time > 22100 && time < 22900) {
              LIBS.translateZ(VIEW_MATRIX, -dt * 0.01);
              LIBS.translateX(VIEW_MATRIX, dt * 0.01);
          }
  
  
  
  
          if (time > 26900 && time < 27700) {
              LIBS.translateZ(VIEW_MATRIX, dt * 0.01);
              LIBS.translateX(VIEW_MATRIX, -dt * 0.01);
          }
  
  
  
  
          // JAKE
          if (time > 24100 && time < 24800) {
              JAKE_MOVE();
              LIBS.rotateY(MODEL_MATRIX_JAKE, -dt * 0.001);
          } else if (time > 25200 && time < 25800) {
              JAKE_MOVE();
              LIBS.rotateY(MODEL_MATRIX_JAKE, dt * 0.001);
          } else if (time > 25900 && time < 41700) {
              JAKE_MOVE();
              LIBS.translateZ(MODEL_MATRIX_JAKE, -dt * 0.001);
          } else if (time > 41700 && time < 43000) {
              JAKE_MOVE();
              LIBS.rotateY(MODEL_MATRIX_JAKE, -dt * 0.001);
          } else if (time > 43000 && time < 58500) {
              JAKE_MOVE();
              LIBS.translateX(MODEL_MATRIX_JAKE, dt * 0.001);
          } else if (time > 58500 && time < 59600) {
              JAKE_MOVE();
              LIBS.translateX(MODEL_MATRIX_JAKE, dt * 0.001);
          } else if (time > 59600 && time < 62800) {
              LIBS.rotateX(MODEL_LEFT_ARM_MATRIX_JAKE, -dt * 0.001);
  
  
  
  
  
  
  
  
          } else if (time > 62800 && time < 64130) {
              LIBS.scaleHeight(MODEL_LEFT_ARM_MATRIX_JAKE, 1.01947);
              if (time > 63600 && time < 63660) {
                  LIBS.rotateX(MODEL_LEFT_ARM_MATRIX_JAKE, dt * 0.01264);
              }
          } if (time > 64100 && time < 73300) {
              if (time > 64100 && time < 64800) {
                  LIBS.translateX(MODEL_KALENG2_MATRIX, -dt * 0.001);
              }
              LIBS.translateZ(MODEL_KALENG2_MATRIX, dt * 0.001);
          } if (time > 65000 && time < 66330) {
              LIBS.scaleHeight(MODEL_LEFT_ARM_MATRIX_JAKE, 1 / 1.01947);
              if (time > 65800 && time < 65860) {
                  LIBS.rotateX(MODEL_LEFT_ARM_MATRIX_JAKE, -dt * 0.01264);
              }
          } if (time > 66330 && time < 69300) {
              LIBS.rotateX(MODEL_LEFT_ARM_MATRIX_JAKE, dt * 0.001);
          } else if (time > 69300 && time < 70350) {
              JAKE_MOVE();
              LIBS.rotateY(MODEL_MATRIX_JAKE, dt * 0.001)
              if (time > 69300 && time < 69910) {
                  LIBS.translateY(MODEL_LEFT_ARM_MATRIX_JAKE, dt * 0.001)
              }
          }
  
  
  
  
          // GUNTER
          if (time > 77600 && time < 88400) {
              LIBS.translateX(VIEW_MATRIX, -dt * 0.001);
          } else if (time > 88400 && time < 95300) {
              GUNTER_MOVE();
              LIBS.translateX(GUNTER_MODEL_MATRIX, -dt * 0.001);
          } else if (time > 95300 && time < 104750) {
              GUNTER_MOVE();
              LIBS.translateY(VIEW_MATRIX, -dt * 0.001);
              LIBS.translateX(GUNTER_MODEL_MATRIX, -dt * 0.00017);
              LIBS.translateY(GUNTER_MODEL_MATRIX, dt * 0.001);
          } else if (time > 105300 && time < 107200) {
              GUNTER_MOVE();
              LIBS.translateX(GUNTER_MODEL_MATRIX, -dt * 0.00082);
          } else if (time > 107200 && time < 108500) {
              GUNTER_MOVE();
              LIBS.rotateY(GUNTER_MODEL_MATRIX, dt * 0.001);
          } else if (time > 109600 && time < 111500) {
              GUNTER_MOVE();
              LIBS.translateZ(GUNTER_MODEL_MATRIX, dt * 0.001);
          } else if (time > 111800 && time < 112000) {
              LIBS.rotateZ(GUNTER_HAND_MODEL_MATRIX2, -dt * 0.001);
          } if (time > 112000 && time < 121200) {
              if (time > 112000 && time < 112800) {
                  LIBS.translateX(MODEL_KALENG1_MATRIX, -dt * 0.001);
                  LIBS.translateX(VIEW_MATRIX, dt * 0.001);
              }
              LIBS.translateZ(MODEL_KALENG1_MATRIX, dt * 0.001);
              LIBS.translateY(VIEW_MATRIX, dt * 0.001);
          } if (time > 121200 && time < 130400) {
              if (time > 121200 && time < 122000) {
                  LIBS.translateX(VIEW_MATRIX, -dt * 0.001);
              }
              LIBS.translateY(VIEW_MATRIX, -dt * 0.001);
          } else if (time > 130400 && time < 133500) {
              GUNTER_MOVE();
              LIBS.rotateY(GUNTER_MODEL_MATRIX, dt * 0.001);
          } else if (time > 133500 && time < 135400) {
              GUNTER_MOVE();
              LIBS.translateZ(GUNTER_MODEL_MATRIX, -dt * 0.001);
          } else if (time > 135400 && time < 137100) {
              GUNTER_MOVE();
              LIBS.rotateY(GUNTER_MODEL_MATRIX, dt * 0.001);
  
  
  
  
          } else if (time > 137100 && time < 139000) {
              GUNTER_MOVE();
              LIBS.translateX(GUNTER_MODEL_MATRIX, dt * 0.00082);
  
  
  
  
          } else if (time > 139000 && time < 148350) {
              GUNTER_MOVE();
              LIBS.translateY(VIEW_MATRIX, dt * 0.001);
              LIBS.translateX(GUNTER_MODEL_MATRIX, dt * 0.00017);
              LIBS.translateY(GUNTER_MODEL_MATRIX, -dt * 0.001);
          }
        
        var axisStart = [0, 0, 0];
        var axisEnd = [1, 1, 1];
        // Menentukan sudut rotasi dalam radian
        var rotationAngle = 270;
        LIBS.rotateAboutAxis(MODEL_KINCIR1_MATRIX, rotationAngle, [120,-40,0], axisStart);
        LIBS.rotateAboutAxis(MODEL_KINCIR2_MATRIX, rotationAngle, [2,-1,0], axisStart);
        LIBS.rotateAboutAxis(MODEL_KINCIR3_MATRIX, rotationAngle, [5,0,1], axisStart);

        if (time > 149350 && time < 153500) {
            LIBS.translateZ(VIEW_MATRIX, -dt * 0.01);
        }
        if (time > 153500 && time < 154500) {
            LIBS.translateY(VIEW_MATRIX, -dt * 0.01);
        }
        if (time > 154500) {
            LIBS.rotateY(VIEW_MATRIX, dt * 0.00080);
        }
  
  
        var radius = 70;
        var pos_x = radius * Math.cos(time * 0.0005);
        var pos_y = radius * Math.sin(time * 0.0005);
        LIBS.set_position(MODEL_MATAHARI_MATRIX, pos_x, 20, pos_y);
  
  
        requestAnimationFrame(animate);
  
  
        GL.flush();
    };
  
  
    animate(0);
  };
  window.addEventListener('load', main);
  
  

