



var innovationNumber = 50;

// Retrieve the canvas information
var canvas = document.getElementById('worldCanvas');
var ctx = canvas.getContext('2d');


var ground = createGround(0, 40, 10000, 40);

var orgs = []

//setup debug draw
var debugDraw = new b2DebugDraw();
debugDraw.SetSprite(ctx);
debugDraw.SetDrawScale(WORLD_DRAW_SCALE);
debugDraw.SetFillAlpha(0.3);
debugDraw.SetLineThickness(1.0);
debugDraw.SetFlags(b2DebugDraw.e_shapeBit | b2DebugDraw.e_jointBit);
world.SetDebugDraw(debugDraw);//*/


var batchSize = 40;

var pool = new GenomePool(batchSize, function() {
	var genes = [];
	for (var i=0; i<batchSize;i++) {
		orgs[i] = CreateBoxOrganism(null);
		genes.push(orgs[i].getGene());
	}
})

var totalMax = 0;
var gen = 0;
var updating = false;
var resetCounter = 30000;
/*var restart = function() {

	var max = -1;
	for (var i=0; i<batchSize;i++) {
		pool.assignFitness(orgs[i].getGene(), orgs[i].getFitness());
		max = Math.max(orgs[i].getFitness(), max);
		//orgs[i] = null;
	}
	if (totalMax < max) {
		totalMax = max;
	}
	console.log("Gen: "+ (gen++) +" | Current best: " + max + " | All time best: " + totalMax + " | Total species: " + pool.getTotalSpecies());

	var newGenes = pool.endTesting();
	if (newGenes != null) {
		for (var i=0; i<batchSize;i++) {
			orgs[i] = CreateBoxOrganism(newGenes[i]);
		}
	}

	pool.startTesting();
	newGenes = null;
};*/

//

/*setInterval(function() {
	//console.log(resetCounter);
	if (resetCounter <= 0) {
		mustReset = true;
	}
	resetCounter -= 1000;
}, 1000);*/

var update = function() {

    var alive = 0;

    for (var i=0; i<batchSize;i++) {
    	if (orgs[i].isAlive()) {
    		orgs[i].update();
    		alive += 1;
    	}
    }

    if (alive == 0) {
		 var max = -1;
		 var newGenes = [];
		for (var i=0; i<batchSize;i++) {
			orgs[i].kill();
			pool.assignFitness(orgs[i].getGene(), orgs[i].getFitness());
			max = Math.max(orgs[i].getFitness(), max);
			//newGenes.push(orgs[i].getGene());
		}
		if (totalMax < max) {
			totalMax = max;
		}
		console.log("Gen: "+ (gen++) +" | Current best: " + max + " | All time best: " + totalMax + " | Total species: " + pool.getTotalSpecies());

		newGenes = pool.endTesting();
		if (newGenes != null) {
			for (var i=0; i<batchSize;i++) {
				orgs[i] = CreateBoxOrganism(newGenes[i].clone());
			}
		}

		pool.startTesting();
		newGenes = null;
    }

    //ctx.save();
	//ctx.clearRect(0, 0, canvas.width, canvas.height);
    //ctx.translate(-200, 0);
    world.Step(
  			(1 / 60)  //frame-rate
        ,  20       //velocity iterations
        ,  20       //position iterations
    );
    world.DrawDebugData();
    world.ClearForces();
    //ctx.restore();

    setTimeout(update, 1000 * (1 / 60));
};

update();




