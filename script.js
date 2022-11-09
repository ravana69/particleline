/**
 * 跨浏览器RAF
 */
window.requestAnimationFrame = (function () {
  return window.requestAnimationFrame ||
    window.webkitRequestAnimationFrame ||
    window.mozRequestAnimationFrame ||
    window.oRequestAnimationFrame ||
    window.msRequestAnimationFrame ||
    function (callback) {
      window.setTimeout(callback, 1000 / 60);
    };
})();

/**
 * 跨浏览器获取3D绘图上下文
 */
let getWebGLContext = function(canvas, opt) {
  let namesArr = ["webgl", "experimental-webgl", "webkit-3d", "moz-webgl"];
  let context = null;
  for (let i = 0; i < namesArr.length; i++) {
    try {
      context = canvas.getContext(namesArr[i], opt);
    } catch(e) {}
    if (context) {
      break;
    }
  }
  return context;
}

function setViewPort(canvas, gl) {
    let b = document.body;
    let d = document.documentElement;
    let fullw = Math.max(b.clientWidth, b.scrollWidth, d.scrollWidth, d.clientWidth);
    let fullh = Math.max(b.clientHeight, b.scrollHeight, d.scrollHeight, d.clientHeight);
    canvas.width = fullw;
    canvas.height = fullh;
    gl.clearColor(1.0, 1.0, 1.0, 1.0);
    gl.viewport(0, 0, canvas.width, canvas.height);
}

/**
 * 编译GLSL ES代码，创建和初始化着色器（顶点和片元）
 * @param {WebGLContext} gl WebGL绘制上下文
 * @param {String} vertexShader 顶点着色器源码
 * @param {String} fragmentShader 片元着色器源码
 * @return program 返回所创建的程序对象
 */
function initShaders(gl, vertexShader, fragmentShader) {
    let vshader;
    vshader = gl.createShader(gl.VERTEX_SHADER);//创建着色器对象
    gl.shaderSource(vshader, vertexShader);//向着色器对象中填充着色器程序的源代码
    gl.compileShader(vshader);//编译着色器
    if (!gl.getShaderParameter(vshader, gl.COMPILE_STATUS)) {//检查着色器的状态
        let msg = "Couldn't compile the vertex shader. The error log is:"
            + "<pre>" + gl.getShaderInfoLog(vshader) + "</pre>";
        //弹出报错信息
        alert(msg);
        gl.deleteShader(vshader);//删除（不是马上删除，会等到程序对象不再使用该着色器后再删除）着色器对象
    return;
        return -1;
    }

    //片元着色器也一样
    let fshader;
    fshader = gl.createShader(gl.FRAGMENT_SHADER);
    gl.shaderSource(fshader, fragmentShader);
    gl.compileShader(fshader);
    if (!gl.getShaderParameter(fshader, gl.COMPILE_STATUS)) {
        let msg = "Couldn't compile the fragment shader.  The error log is:"
            + "<pre>" + gl.getShaderInfoLog(fshader) + "</pre>";
        alert(msg);
        gl.deleteShader(fshader);
        return -1;
    }

    //创建程序对象
    let program = gl.createProgram();
    //为程序对象分配着色器
    gl.attachShader(program, vshader);
    gl.attachShader(program, fshader);
    //连接程序对象(在为程序对象分配了两个着色器对象后，还需要将顶点着色器和片元着色器连接起来)
    gl.linkProgram(program);

    if (!gl.getProgramParameter(program, gl.LINK_STATUS)) {
        let msg = "Unable to initialise shaders.  The error log is:"
            + "<pre>" + gl.getProgramInfoLog(program) + "</pre>";
        alert(msg);
        gl.deleteProgram(gl.program);
        gl.deleteProgram(vshader);
        gl.deleteProgram(fshader);
        return -1;
    }

    return program;
}

class LGMV {
    _argumentsToArray(args) {
        return [].concat.apply([], Array.from(args));
    }
    radians( degrees ) {
        return degrees * Math.PI / 180.0;
    }
    //向量初始化
    vec2() {
        let result = this._argumentsToArray(arguments);
        switch (result.length) {
            case 0: result.push(0.0);
            case 1: result.push(0.0);
        }
        return result.splice(0, 2);
    }
    vec3() {
        let result = this._argumentsToArray(arguments);
        switch (result.length) {
            case 0: result.push(0.0);
            case 1: result.push(0.0);
            case 2: result.push(0.0);
        }
        return result.splice(0, 3);
    }
    vec4() {
        let result = this._argumentsToArray(arguments);
        switch (result.length) {
            case 0: result.push(0.0);
            case 1: result.push(0.0);
            case 2: result.push(0.0);
            case 3: result.push(1.0);
        }
        return result.splice(0, 4);
    }
    //矩阵初始化
    mat2() {
        let arg = this._argumentsToArray(arguments);
        let m = [];
        switch (arg.length) {
            case 0:
                arg[0] = 1;
            case 1:
                m = [
                    this.vec2(arg[0], 0.0),
                    this.vec2(0.0, arg[0])
                ];
                break;
            default:
                m.push(this.vec2(arg)); arg.splice(0, 2);
                m.push(this.vec2(arg));
                break;
        }
        m.matrix = true;
        return m;
    }
    mat3() {
        let arg = this._argumentsToArray(arguments);
        let m = [];
        switch (arg.length) {
            case 0:
                arg[0] = 1;
            case 1:
                m = [
                    this.vec3(arg[0], 0.0, 0.0),
                    this.vec3(0.0, arg[0], 0.0),
                    this.vec3(0.0, 0.0, arg[0])
                ];
                break;
            default:
                m.push(this.vec3(arg)); arg.splice(0, 3);
                m.push(this.vec3(arg)); arg.splice(0, 3);
                m.push(this.vec3(arg));
                break;
        }
        m.matrix = true;
        return m;
    }
    mat4() {
        let arg = this._argumentsToArray(arguments);
        let m = [];
        switch (arg.length) {
            case 0:
                arg[0] = 1;
            case 1:
                m = [
                    this.vec4(arg[0], 0.0, 0.0, 0.0),
                    this.vec4(0.0, arg[0], 0.0, 0.0),
                    this.vec4(0.0, 0.0, arg[0], 0.0),
                    this.vec4(0.0, 0.0, 0.0, arg[0])
                ];
                break;

            default:
                m.push(this.vec4(arg)); arg.splice(0, 4);
                m.push(this.vec4(arg)); arg.splice(0, 4);
                m.push(this.vec4(arg)); arg.splice(0, 4);
                m.push(this.vec4(arg));
                break;
        }
        m.matrix = true;
        return m;
    }
    //顶点方法
    //点乘
    dot(u, v) {
        if (u.length != v.length) {
            throw "this.dot(): vectors are not the same dimension";
        }

        let sum = 0.0;
        for (let i = 0; i < u.length; ++i) {
            sum += u[i] * v[i];
        }

        return sum;
    }
    //叉乘
    cross(u, v) {
        if (!Array.isArray(u) || u.length < 3) {
            throw "this.cross(): first argument is not a vector of at least 3";
        }

        if (!Array.isArray(v) || v.length < 3) {
            throw "this.cross(): second argument is not a vector of at least 3";
        }

        let result = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0]
        ];

        return result;
    }
    negate(u) {
        let result = [];
        for (let i = 0; i < u.length; ++i) {
            result.push(-u[i]);
        }

        return result;
    }
    //矩阵方法
    //相加
    add(u, v) {
        let result = [];

        if (u.matrix && v.matrix) {
            if (u.length != v.length) {
                throw "add(): trying to add matrices of different dimensions";
            }

            for (let i = 0; i < u.length; ++i) {
                if (u[i].length != v[i].length) {
                    throw "add(): trying to add matrices of different dimensions";
                }
                result.push([]);
                for (let j = 0; j < u[i].length; ++j) {
                    result[i].push(u[i][j] + v[i][j]);
                }
            }

            result.matrix = true;

            return result;
        }
        else if (u.matrix && !v.matrix || !u.matrix && v.matrix) {
            throw "add(): trying to add matrix and non-matrix variables";
        }
        else {
            if (u.length != v.length) {
                throw "add(): vectors are not the same dimension";
            }

            for (let i = 0; i < u.length; ++i) {
                result.push(u[i] + v[i]);
            }

            return result;
        }
    }
    //相减
    subtract(u, v) {
        let result = [];

        if (u.matrix && v.matrix) {
            if (u.length != v.length) {
                throw "subtract(): trying to subtract matrices" +
                " of different dimensions";
            }

            for (let i = 0; i < u.length; ++i) {
                if (u[i].length != v[i].length) {
                    throw "subtract(): trying to subtact matrices" +
                    " of different dimensions";
                }
                result.push([]);
                for (let j = 0; j < u[i].length; ++j) {
                    result[i].push(u[i][j] - v[i][j]);
                }
            }

            result.matrix = true;

            return result;
        }
        else if (u.matrix && !v.matrix || !u.matrix && v.matrix) {
            throw "subtact(): trying to subtact  matrix and non-matrix variables";
        }
        else {
            if (u.length != v.length) {
                throw "subtract(): vectors are not the same length";
            }

            for (let i = 0; i < u.length; ++i) {
                result.push(u[i] - v[i]);
            }

            return result;
        }
    }
    //相乘
    mult(u, v) {
        let result = [];

        if (u.matrix && v.matrix) {
            if (u.length != v.length) {
                throw "mult(): trying to add matrices of different dimensions";
            }

            for (let i = 0; i < u.length; ++i) {
                if (u[i].length != v[i].length) {
                    throw "mult(): trying to add matrices of different dimensions";
                }
            }

            for (let i = 0; i < u.length; ++i) {
                result.push([]);

                for (let j = 0; j < v.length; ++j) {
                    let sum = 0.0;
                    for (let k = 0; k < u.length; ++k) {
                        sum += u[i][k] * v[k][j];
                    }
                    result[i].push(sum);
                }
            }

            result.matrix = true;

            return result;
        }

        if (u.matrix && (u.length == v.length)) {
            for (let i = 0; i < v.length; i++) {
                let sum = 0.0;
                for (let j = 0; j < v.length; j++) {
                    sum += u[i][j] * v[j];
                }
                result.push(sum);
            }
            return result;
        }



        else {
            if (u.length != v.length) {
                throw "mult(): vectors are not the same dimension";
            }

            for (let i = 0; i < u.length; ++i) {
                result.push(u[i] * v[i]);
            }

            return result;
        }
    }
    //相等
    equal(u, v) {
        if (u.length != v.length) { return false; }

        if (u.matrix && v.matrix) {
            for (let i = 0; i < u.length; ++i) {
                if (u[i].length != v[i].length) { return false; }
                for (let j = 0; j < u[i].length; ++j) {
                    if (u[i][j] !== v[i][j]) { return false; }
                }
            }
        }
        else if (u.matrix && !v.matrix || !u.matrix && v.matrix) {
            return false;
        }
        else {
            for (let i = 0; i < u.length; ++i) {
                if (u[i] !== v[i]) { return false; }
            }
        }

        return true;
    }
    //归一化
    normalize(u, excludeLastComponent) {
        if (excludeLastComponent) {
            let last = u.pop();
        }

        let len = this.length(u);

        if (!isFinite(len)) {
            throw "normalize: vector " + u + " has zero length";
        }

        for (let i = 0; i < u.length; ++i) {
            u[i] /= len;
        }

        if (excludeLastComponent) {
            u.push(last);
        }

        return u;
    }
    transpose(m) {
        if (!m.matrix) {
            return "transpose(): trying to transpose a non-matrix";
        }

        let result = [];
        for (let i = 0; i < m.length; ++i) {
            result.push([]);
            for (let j = 0; j < m[i].length; ++j) {
                result[i].push(m[j][i]);
            }
        }

        result.matrix = true;

        return result;
    }
    //长度
    length(u) {
        return Math.sqrt(this.dot(u, u));
    }
    //数组类型转换
    flatten(v) {
        if (v.matrix === true) {
            //矩阵先转置（WebGL使用列主序）
            v = this.transpose(v);
        }
        let n = v.length;
        let elemsAreArrays = false;

        if (Array.isArray(v[0])) {
            elemsAreArrays = true;
            n *= v[0].length;
        }

        let floats = new Float32Array(n);

        if (elemsAreArrays) {
            let idx = 0;
            for (let i = 0; i < v.length; ++i) {
                for (let j = 0; j < v[i].length; ++j) {
                    floats[idx++] = v[i][j];
                }
            }
        }
        else {
            for (let i = 0; i < v.length; ++i) {
                floats[i] = v[i];
            }
        }

        return floats;
    }
    mix(u, v, s) {
        if (typeof s !== "number") {
            throw "mix: the last paramter " + s + " must be a number";
        }

        if (u.length != v.length) {
            throw "vector dimension mismatch";
        }

        let result = [];
        for (let i = 0; i < u.length; ++i) {
            result.push((1.0 - s) * u[i] + s * v[i]);
        }

        return result;
    }
    //MVP相关
    lookAt(eye, at, up) {
        if (!Array.isArray(eye) || eye.length != 3) {
            throw "lookAt(): first parameter [eye] must be an a this.vec3";
        }

        if (!Array.isArray(at) || at.length != 3) {
            throw "lookAt(): first parameter [at] must be an a this.vec3";
        }

        if (!Array.isArray(up) || up.length != 3) {
            throw "lookAt(): first parameter [up] must be an a this.vec3";
        }

        if (this.equal(eye, at)) {
            return this.this.mat4();
        }

        let v = this.normalize(this.subtract(at, eye));  // view direction vector
        let n = this.normalize(this.cross(v, up));       // perpendicular vector
        let u = this.normalize(this.cross(n, v));        // "new" up vector

        v = this.negate(v);

        let result = this.mat4(
            this.vec4(n, -this.dot(n, eye)),
            this.vec4(u, -this.dot(u, eye)),
            this.vec4(v, -this.dot(v, eye)),
            this.vec4()
        );
        return result;
    }

    //平行投影
    ortho(left, right, bottom, top, near, far) {
        //视见体是对称的
        //坐标映射：x => [-1.1] y => [-1,1] z => [0,1]
        //矩阵表示(通用): 
        //r => right l => left t => top b => bottom f => far n => near 
        //[ 2/(r-l), 0, 0, -(r+l)/(r-l)
        //  0, 2/(t-b), 0, -(t+b)/(t-b)
        //  0, 0, -2/(f-n), -(f+n)/(f-n)
        //  0, 0, 0, 1]
        // 对称时：r = -l, t = -b ....
        // 可简化成：
        //[ 1/r, 0, 0, 0
        //  0, 1/t, 0, 0
        //  0, 0, -2/(f-n), -(f+n)/(f-n)
        //  0, 0, 0, 1]
         if (left == right) { throw "ortho(): left and right are equal"; }
        if (bottom == top) { throw "ortho(): bottom and top are equal"; }
        if (near == far) { throw "ortho(): near and far are equal"; }

        let w = right - left;
        let h = top - bottom;
        let d = far - near;

        let result = this.mat4();
        result[0][0] = 2.0 / w;
        result[1][1] = 2.0 / h;
        result[2][2] = -2.0 / d;
        result[0][3] = -(left + right) / w;
        result[1][3] = -(top + bottom) / h;
        result[2][3] = -(near + far) / d;

        return result;
    }

    //透视投影
    perspective(fovy, aspect, near, far) {
        //aspect = width /height => 2r/2t(对称的情况下) => r/t
        //f = 1/ (top/near)  视角还少一半 => n/t
        var f = 1.0 / Math.tan(this.radians(fovy) / 2);
        var d = far - near;
        //和平息投影原理时一样的，后面的推导加个相似三角形
        //通用：
        //[ 2n/(r-l), 0, (r+l)/(r-l), 0
        //  0, 2n/(t-b), (t+b)/(t-b), 0
        //  0, 0, -(f+n)/(f-n), -2fn/(f-n)
        //  0, 0, -1, 0]
        //对称简化后：
        //[ n/r, 0, 0, 0
        //  0, n/t, 0, 0
        //  0, 0, -(f+n)/(f-n), -2fn/(f-n)
        //  0, 0, -1, 0 ]
        var result = this. mat4();
        result[0][0] = f / aspect;
        result[1][1] = f;
        result[2][2] = -(near + far) / d;//这里的符号问题
        result[2][3] = -2 * near * far / d;//这里的符号问题
        result[3][2] = -1;
        result[3][3] = 0.0;//默认是1，手动置0

        return result;
    }
}

//起手式
let canvas = document.getElementById("webgl");
let gl = getWebGLContext(canvas);
setViewPort(canvas, gl);
let MV = new LGMV();

let ratio,
    vertices,//顶点坐标
    velocities,
    freqArr,
    cw,
    ch,
    colorLoc,
    thetaArr,
    velThetaArr,
    velRadArr,
    boldRateArr,
    drawType,
    numLines = 40000;
var target = [];
var randomTargetXArr = [], randomTargetYArr = [];
drawType = 2;

gl.clearColor(0.0, 0.0, 0.0, 1.0);
gl.clearDepth(1.0);
gl.enable(gl.BLEND);//开启gl.BLEND(混合)
//gl.disable(gl.DEPTH_TEST);//关闭gl.DEPTH_TEST(隐藏面消除)
gl.blendFunc(gl.SRC_ALPHA, gl.ONE);//混合的模式
gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);

setup();

function setup() {
    vertices = [];
    velThetaArr = [];
    velRadArr = [];
    ratio = cw / ch;
    velocities = [];
    thetaArr = [];
    freqArr = [];
    boldRateArr = [];

    for (var ii = 0; ii < numLines; ii++) {
        var rad = (0.1 + .2 * Math.random());
        var theta = Math.random() * Math.PI * 2;
        var velTheta = Math.random() * Math.PI * 2 / 30;
        var freq = Math.random() * 0.12 + 0.03;
        var boldRate = Math.random() * .04 + .01;
        var randomPosX = (Math.random() * 2 - 1) * window.innerWidth / window.innerHeight;
        var randomPosY = Math.random() * 2 - 1;

        vertices.push(rad * Math.cos(theta), rad * Math.sin(theta), 1.83);
        vertices.push(rad * Math.cos(theta), rad * Math.sin(theta), 1.83);

        thetaArr.push(theta);
        velThetaArr.push(velTheta);
        velRadArr.push(rad);
        freqArr.push(freq);
        boldRateArr.push(boldRate);

        randomTargetXArr.push(randomPosX);
        randomTargetYArr.push(randomPosY);
    }

    freqArr = new Float32Array(freqArr);
    vertices = new Float32Array(vertices);
    velocities = new Float32Array(velocities);

    thetaArr = new Float32Array(thetaArr);
    velThetaArr = new Float32Array(velThetaArr);
    velRadArr = new Float32Array(velRadArr);
}

let vshader = document.getElementById("vertex-shader").text;
let fshader = document.getElementById("fragment-shader").text;
let program = initShaders(gl, vshader, fshader);
gl.useProgram(program);

let vBuffer = gl.createBuffer();
gl.bindBuffer(gl.ARRAY_BUFFER, vBuffer);
gl.bufferData(gl.ARRAY_BUFFER, vertices, gl.STATIC_DRAW);

let vPosition = gl.getAttribLocation(program, "vPosition");
gl.vertexAttribPointer(vPosition, 3, gl.FLOAT, false, 0, 0);
gl.enableVertexAttribArray(vPosition);

let modelView = gl.getUniformLocation(program, "modelView");
let projection = gl.getUniformLocation(program, "projection");

//透视投影
var fieldOfView = 30.0;
var aspectRatio = canvas.width / canvas.height;
var nearPlane = 1.0;
var farPlane = 10000.0;
var pMatrix = MV.perspective(fieldOfView , aspectRatio, nearPlane, farPlane);
gl.uniformMatrix4fv(projection, false, MV.flatten(pMatrix));

const at = MV.vec3(0.0, 0.0, 0.0);//look at point
const up = MV.vec3(0.0, 1.0, 0.0);//视点上方向
eye = MV.vec3(0, 0, -2);
    
var MVMatrix = MV.lookAt(eye, at, up);
gl.uniformMatrix4fv(modelView, false, MV.flatten(MVMatrix));

animate();//这是核心
setTimeout(timer, 1500);

var count = 0;
var cn = 0;

function animate() {
    requestAnimationFrame(animate);
    drawScene();
}

function drawScene() {
    draw();
    gl.lineWidth(1);
    gl.bufferData(gl.ARRAY_BUFFER, vertices, gl.DYNAMIC_DRAW);
    gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
    gl.drawArrays(gl.LINES, 0, numLines);

    gl.flush();
}

function draw() {
    switch (drawType) {
        case 0:
            draw1();
            break;
        case 1:
            draw2();
            break;
        case 2:
            draw0();
            break;
    }
}

function draw0() {

    var i, n = vertices.length, p, bp;
    var px, py;
    var pTheta;
    var rad;
    var num;
    var targetX, targetY;

    for (i = 0; i < numLines * 2; i += 2) {//因为是线，粒子数量乘以2，然后将多复制出来的一份粒子坐标略做偏移，之后相连，效果就出来了？？？？？并不对i+=2才是导致numLines * 2的原因
        count += .3;
        bp = i * 3;

        vertices[bp] = vertices[bp + 3];
        vertices[bp + 1] = vertices[bp + 4];//确保线起始点一致

        num = parseInt(i / 2);
        targetX = randomTargetXArr[num];//init时的随机位置
        targetY = randomTargetYArr[num];


        px = vertices[bp + 3];
        px += (targetX - px) * (Math.random() * .04 + .6);//相较于上一时刻的位置略做偏移
        vertices[bp + 3] = px;//所以虽然是线但是两个点的坐标是完全一致的，除了亮度增加以外看着就是一个点了


        //py = (Math.sin(cn) + 1) * .2 * (Math.random() * .5 - .25);
        py = vertices[bp + 4];
        py += (targetY - py) * (Math.random() * .04 + .06);//位置略做偏移
        vertices[bp + 4] = py;

    }
}

function draw1() {

    var i, n = vertices.length, p, bp;
    var px, py;
    var pTheta;
    var rad;
    var num;
    var targetX, targetY;

    for (i = 0; i < numLines * 2; i += 2) {
        count += .3;
        bp = i * 3;

        vertices[bp] = vertices[bp + 3];
        vertices[bp + 1] = vertices[bp + 4];

        num = parseInt(i / 2);
        pTheta = thetaArr[num];
        rad = velRadArr[num];

        pTheta += velThetaArr[num];
        thetaArr[num] = pTheta;

        targetX = rad * 2 * Math.cos(pTheta);
        targetY = rad * 2 * Math.sin(pTheta);

        px = vertices[bp + 3];
        px += (targetX - px) * (Math.random() * .1 + .1);//朝方向缓动
        vertices[bp + 3] = px;

        //py = (Math.sin(cn) + 1) * .2 * (Math.random() * .5 - .25);
        py = vertices[bp + 4];
        py += (targetY - py) * (Math.random() * .1 + .1);
        vertices[bp + 4] = py;
    }
}

function draw2() {
    cn += .1;

    var i, n = vertices.length, p, bp;
    var px, py;
    var pTheta;
    var rad;
    var num;

    for (i = 0; i < numLines * 2; i += 2) {
        count += .3;
        bp = i * 3;
        // copy old positions

        vertices[bp] = vertices[bp + 3];
        vertices[bp + 1] = vertices[bp + 4];

        num = parseInt(i / 2);
        pTheta = thetaArr[num];

        rad = velRadArr[num];// + Math.cos(pTheta + i * freqArr[i]) *  boldRateArr[num];

        pTheta += velThetaArr[num];
        thetaArr[num] = pTheta;

        px = vertices[bp + 3];
        px += rad * Math.cos(pTheta) * 0.5;//cos、sin都会有负值
        vertices[bp + 3] = px;


        //py = (Math.sin(cn) + 1) * .2 * (Math.random() * .5 - .25);
        py = vertices[bp + 4];

        py += rad * Math.sin(pTheta) * 0.5;
        //p *= ( Math.random() -.5);
        vertices[bp + 4] = py;
    }
}

function timer() {
    drawType = (drawType + 1) % 3;

    setTimeout(timer, 1500);
}