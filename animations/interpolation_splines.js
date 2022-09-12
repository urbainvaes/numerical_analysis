let xmin = -1;
let xmax = 1;
let ymax = 1.2;
let ymin = -.5;
let board = JXG.JSXGraph.initBoard('mybox', {boundingbox:[xmin,ymax,xmax,ymin], axis:false,
    showCopyright:false, grid:false, showFullscreen:false, ShowNavigation:false, pan: {enabled:false}, zoom:{enabled:true}});

let p = [];
let runge = function(x) {return 1 / (1 + 25*x*x);}
let cos = function(x) {return (Math.cos(2*Math.PI*x) + 1)/2;}
let currentf = cos;
let ipoint = 0;
let npoints = 6;
for (; ipoint<npoints; ipoint++) {
    x = xmin + (xmax-xmin) * (ipoint+1) / ((npoints+2)-1);
    y = currentf(x);
    p[ipoint] = board.create('point', [x,y], {size: 4, face: 'o', withLabel:false, showInfobox:false, highlightFillColor:'black', highlightStrokeColor:'black'});
}

let f = JXG.Math.Numerics.lagrangePolynomial(p);
var f_txt = JXG.Math.Numerics.lagrangePolynomialTerm(p, 2, 'x', ' * ');
let fspline = JXG.Math.Numerics.CardinalSpline(p, .5, "uniform");
let interp = board.create('functiongraph', [f,xmin, xmax], {highlight:false, strokeColor:'blue', strokeOpacity:0.5, strokeWidth:5, label:"Lagrange interpolant"});
let interp_spline = board.create('spline', p, {highlight:false, strokeColor:'green', strokeOpacity:0.3, strokeWidth:5});

let getMouseCoords = function(e, i) {
    var cPos = board.getCoordsTopLeftCorner(e, i),
        absPos = JXG.getPosition(e, i),
        dx = absPos[0]-cPos[0],
        dy = absPos[1]-cPos[1];

    return new JXG.Coords(JXG.COORDS_BY_SCREEN, [dx, dy], board);
},
    down = function(e) {
        var canCreate = true, i, coords, el;

        if (e[JXG.touchProperty]) {
            // index of the finger that is used to extract the coordinates
            i = 0;
        }
        coords = getMouseCoords(e, i);

        for (el in board.objects) {
            if(JXG.isPoint(board.objects[el]) && board.objects[el].hasPoint(coords.scrCoords[1], coords.scrCoords[2])) {
                canCreate = false;
                break;
            }
        }

        if (canCreate) {
            x = coords.usrCoords[1];
            y = coords.usrCoords[2];
            p[ipoint] = board.create('point', [x, y], {size: 4, face: 'o', withLabel:false, showInfobox:false, highlightFillColor:'black', highlightStrokeColor:'black'});
            ipoint = ipoint+1;
        }
        board.update()
    }
board.on('down', down);
