
var DEGTORAD = Math.PI / 180.0;

// Creaate shortcut names for commons
var b2Vec2 = Box2D.Common.Math.b2Vec2,
	b2BodyDef = Box2D.Dynamics.b2BodyDef,
	b2Body = Box2D.Dynamics.b2Body,
	b2FixtureDef = Box2D.Dynamics.b2FixtureDef,
	b2Fixture = Box2D.Dynamics.b2Fixture,
	b2World = Box2D.Dynamics.b2World,
	b2MassData = Box2D.Collision.Shapes.b2MassData,
	b2PolygonShape = Box2D.Collision.Shapes.b2PolygonShape,
	b2CircleShape = Box2D.Collision.Shapes.b2CircleShape,
	b2DebugDraw = Box2D.Dynamics.b2DebugDraw,
	b2WeldJointDef = Box2D.Dynamics.Joints.b2WeldJointDef,
	b2RevoluteJointDef = Box2D.Dynamics.Joints.b2RevoluteJointDef,
	b2ContactListener = Box2D.Dynamics.b2ContactListener;

// 
var world = new b2World(new b2Vec2(0.0, 5), false);

var CANVAS_HEIGHT = 400;

// 
var WORLD_DRAW_SCALE = 30.0;

// The fixture definition for static ground bodies
var groundFixDef = new b2FixtureDef;
groundFixDef.density = 1.0;
groundFixDef.friction = 0.5;
groundFixDef.restitution = 0.2;
groundFixDef.filter.categoryBits = 0x0002;
groundFixDef.filter.maskBits = 0x0004;

// The fixture definition for oganism bodies
var organismFixDef = new b2FixtureDef;
organismFixDef.density = 0.4;
organismFixDef.friction = 0.5;
organismFixDef.restitution = 0.2;
organismFixDef.filter.categoryBits = 0x0004;
organismFixDef.filter.maskBits = 0x0002;

var sensorFixtureDef = new b2FixtureDef;
sensorFixtureDef.density = 0.0001;
sensorFixtureDef.isSensor = true;
sensorFixtureDef.filter.categoryBits = 0x0004;
sensorFixtureDef.filter.maskBits = 0x0002;

// Create the joints
var joint = new b2RevoluteJointDef;
joint.enableMotor=true;
joint.maxMotorTorque=30.01;
joint.motorSpeed=0;
joint.enableLimit = true;
joint.lowerAngle = -45 * DEGTORAD;
joint.upperAngle =  45 * DEGTORAD;
joint.collideConnected = false;

var weldJointDef = new b2WeldJointDef;
weldJointDef.referenceAngle = 0;
weldJointDef.dampingRatio = 0;
weldJointDef.frequencyHz = 0;


var cl = new b2ContactListener;

cl.BeginContact = function(contact) {
	var fixtureA = contact.GetFixtureA();
    var fixtureB = contact.GetFixtureB();
    if (fixtureA.IsSensor()) {
    	fixtureA.SetUserData(true);
    } else if (fixtureB.IsSensor()) {
    	fixtureB.SetUserData(true);
    }
};

cl.EndContact = function(contact) {
	var fixtureA = contact.GetFixtureA();
    var fixtureB = contact.GetFixtureB();
    if (fixtureA.IsSensor()) {
    	fixtureA.SetUserData(false);
    } else if (fixtureB.IsSensor()) {
    	fixtureB.SetUserData(false);
    }
};

world.SetContactListener(cl);

/** Creates a block of ground at the bottom of the world. */
var createGround = function(x1, h1, x2, h2) {
	var bodyDef = new b2BodyDef;
    bodyDef.type = b2Body.b2_staticBody;
    bodyDef.position.x = 0;
    bodyDef.position.y = 0;

    // Set the shape of the ground
    groundFixDef.shape = new b2PolygonShape;
    groundFixDef.shape.SetAsArray([
    	new b2Vec2(x1 / WORLD_DRAW_SCALE, (CANVAS_HEIGHT - h1) / WORLD_DRAW_SCALE),
    	new b2Vec2(x2 / WORLD_DRAW_SCALE, (CANVAS_HEIGHT - h2) / WORLD_DRAW_SCALE),
    	new b2Vec2(x2 / WORLD_DRAW_SCALE, CANVAS_HEIGHT / WORLD_DRAW_SCALE),
    	new b2Vec2(x1 / WORLD_DRAW_SCALE, CANVAS_HEIGHT / WORLD_DRAW_SCALE),
    	], 4);
 
    // Create and retrieve the body 
    var b = world.CreateBody(bodyDef);
    b.CreateFixture(groundFixDef);
    return b;
};

var innov = null;

/** */
var Organism = function(genome, bodies, head, joints, sensors) {
	var MEMORY_NEURONS = 0;
	var MAX_MOTOR_SPEED = 7;
	var INITIAL_HIDDEN_NEURONS = 5;

	// The initial health of the organism
	var health = 500;

	// The current fitness of the organism
	var fitness = 0;

	// 
	var aliveState = "ALIVE";

	//var leg = [];

	// Set the neural network of the organism
	if (genome == null) {

		if (innov == null) {
			innov = new InnovationNumber(100);
		}

		genome = new Genome(3, joints.length, innov);//MEMORY_NEURONS + joints.length + 2 + sensors.length, MEMORY_NEURONS + joints.length);
		genome.addVariableGene(150, 1, 300); // Gait duration
		/*for (var i=0; i<joints.length; i++) {
			genome.addValueGene(0, -1, 1); // Time staggered
			genome.addValueGene(0, -1, 1); // Direction staggered
		}
		var inputNeurons = genome.getInputIds();
		var outputNeurons = genome.getOutputIds();
		var count = inputNeurons.length + outputNeurons.length;
		for (var i = 0; i < INITIAL_HIDDEN_NEURONS; i++) {
			for (var j = 0; j < inputNeurons.length; j++) {
				genome.addGene(new Gene(count, inputNeurons[j], count, Math.random() * 2 - 1))
			}
			for (var j = 0; j < outputNeurons.length; j++) {
				genome.addGene(new Gene(count, count, inputNeurons[j], count, Math.random() * 2 - 1))
			}
		}*/
	}

	// Set the initial values of the memory neurons
	var memory = [];
	for (var i = 0; i < MEMORY_NEURONS; i++) {
		memory.push(0);
	}

	/** */
	this.isAlive = function() {
		return aliveState == "ALIVE";
	} 

	/** Sets the organism to stop simulating. The destruction of the physic bodies will occur during the physic step. */
	this.flagToKill = function() {
		if (aliveState == "ALIVE") {
			aliveState = "TO_DIE";
		}
	}

	this.kill = function() {
		if (aliveState != "DEAD") {
			aliveState = "DEAD";
			for (var i = 0; i < joints.length; i++) {
				for (var j = 0; j < joints[i].length; i++) {
					world.DestroyJoint(joints[i][j]);
				}
			}
			for (var i = 0; i < bodies.length; i++) {
				world.DestroyBody(bodies[i]);
			}

			bodies = null;
			head = null;
			joints = null;
			memory = null;
			sensors = null;
		}
	}

	this.getGene = function() {
		return genome;
	}

	this.getFitness = function() {
		return fitness;
	}

	var timeMoving = 1;
	var timeTotal = 0;

	var gaitTime = 0;

	var gaitDirections = [];
	for (var i=0; i<joints.length;i++) {
		gaitDirections[i] = 1;
	}

	/** */
	this.update = function() {
		if (aliveState == "ALIVE") {

			var inputValues = [];
			inputValues.push(1);
			var v = (gaitTime / genome.getVariableGene(0)) * 2 * Math.PI;
			inputValues.push(Math.sin(v));
			inputValues.push(Math.cos(v))

			// Pass in the memory output neuronn values from the last propagation
			/*for (var i = 0; i < MEMORY_NEURONS; i++) {
				inputValues.push(memory[i]);
			}

			// Add the value of the bias neuron
			inputValues.push(1);
			inputValues.push(head.GetAngle());

			// Add information about each joint to the input values of the neural network
			for (var i = 0; i < joints.length; i++) {
				inputValues.push(joints[i].GetJointAngle());
			}

			for (var i = 0; i < sensors.length; i++) {
				inputValues.push((sensors[i].GetUserData() == true) ? 1 : 0);
			}*/

			// Run the values through the network
			var outputValues = genome.propagate(inputValues);

			// Store the values of the memory node from the previous iteration
			/*for (var i = 0; i < MEMORY_NEURONS; i++) {
				memory[i] = outputValues.pop();
			}*/

			// Set the speed of the joints of the organisms
			var moving = false;
			for (var i = 0; i < joints.length; i++) {
				var amount = (outputValues[i] * 2.0) - 1.0;

				joints[i].SetMotorSpeed(MAX_MOTOR_SPEED * amount);

				/*if (Math.abs(amount) < 0.1) {
					joints[i].SetMotorSpeed(0); // No operation
				} else {
					joints[i].SetMotorSpeed(MAX_MOTOR_SPEED * amount);
					/*
					if (!(amount > 0 && joints[i].GetJointAngle() > (43.0 * DEGTORAD)) && 
						!(amount < 0 && joints[i].GetJointAngle() < (-43.0 * DEGTORAD))) {
						joints[i].SetMotorSpeed(MAX_MOTOR_SPEED * amount);
						moving = true;
					} else {
						// The joints cannot move any farther so do not move the joints
						joints[i].SetMotorSpeed(0);
					}
				}*/
			}

			gaitTime += 1;
			if (gaitTime >= genome.getVariableGene(0)) {
				//gaitDirection *= -1;
				gaitTime = 0;
			}

			// Start killing the organism if it is not moving
			if (moving == false) {
				health -= 1;
			} else {
				timeMoving += 1;
			}

			timeTotal += 1;

			if (!isNaN(head.GetPosition().x) && head.GetPosition().y < 10.5) {
				fitness = (head.GetPosition().x - 2.0) * (timeMoving / timeTotal);

				if (Math.abs(fitness - head.GetPosition().x) < 0.01) {
            		health -= 1;
         		}
			}

			// Kill the organism if it has no health
			if (health <= 0) {
				this.kill();
			} 
		}
	}
};

var CreateBoxOrganism = function(genome) {

	var bodyDef = new b2BodyDef;
    bodyDef.type = b2Body.b2_dynamicBody;
    bodyDef.position.x = 2;
    bodyDef.position.y = 10;

    // Set the bosy fixture shape
    organismFixDef.shape = new b2PolygonShape;
    organismFixDef.shape.SetAsArray([
    	new b2Vec2(0, 0),
    	new b2Vec2(2, 0),
    	new b2Vec2(2, 0.5),
    	new b2Vec2(0, 0.5)
    	], 4);

    // Create the organism main body
    var centralBody = world.CreateBody(bodyDef);
    centralBody.CreateFixture(organismFixDef);

    // Set the sensor fixture shape
    var sensorDef = new b2BodyDef;
	    sensorDef.type = sensorDef.b2_dynamicBody;
	    sensorFixtureDef.shape = new b2PolygonShape;
	    sensorFixtureDef.shape.SetAsArray([
	    	new b2Vec2(-0.01, 0.9),
	    	new b2Vec2(0.26, 0.9),
	    	new b2Vec2(0.26, 1.01),
	    	new b2Vec2(-0.01, 1.01)
	    	], 4);

	var createLeg = function(x, y, jx) {
		// Set the shape of the leg
		organismFixDef.shape.SetAsArray([
	    	new b2Vec2(0, 0),
	    	new b2Vec2(0.25, 0),
	    	new b2Vec2(0.25, 1),
	    	new b2Vec2(0, 1)
    	], 4);


		// Create the upper part of the leg
		bodyDef.position.x = x;
    	bodyDef.position.y = y;
	    var upperLeg = world.CreateBody(bodyDef);
	    upperLeg.CreateFixture(organismFixDef);

	    // Create the lower part of the leg
	    bodyDef.position.x = x;
	    bodyDef.position.y = y + 1;
	    var lowerLeg = world.CreateBody(bodyDef);
	    lowerLeg.CreateFixture(organismFixDef);

	    // Create the sensor for the leg
	    var s1 = lowerLeg.CreateFixture(sensorFixtureDef);

	    // Create the joint that connects the lower and upper leg
	    joint.bodyA = centralBody;
	  	joint.bodyB = upperLeg;
	  	joint.localAnchorA.Set(jx, 0.25);
	    joint.localAnchorB.Set(0.125,0.125);
	    var upperJoint = world.CreateJoint(joint);

	    joint.bodyA = upperLeg;
	  	joint.bodyB = lowerLeg;
	  	joint.localAnchorA.Set(0.125, 0.875);
	    joint.localAnchorB.Set(0.125, 0.125);
		var lowerJoint = world.CreateJoint(joint);

	    return {
	    	"bodies" : [upperLeg, lowerLeg],
	    	"sensors" : [s1],
	    	"joints" : [upperJoint, lowerJoint]
	    }
	}

	var leftFrontLeg = createLeg(2, 10, 0.125);
	var leftLegBack = createLeg(2, 10, 0.125);

	var rightFrontLeg = createLeg(4, 10, 1.875);
	var rightLegBack = createLeg(4, 10, 1.875);

	var bodies = leftFrontLeg.bodies.concat(leftLegBack.bodies).concat(rightFrontLeg.bodies).concat(rightLegBack.bodies);
	var sensors = leftFrontLeg.sensors.concat(leftLegBack.sensors).concat(rightFrontLeg.sensors).concat(rightLegBack.sensors);
	var joints = leftFrontLeg.joints.concat(leftLegBack.joints).concat(rightFrontLeg.joints.concat(rightLegBack.joints));

	var frontLegJoints = [leftFrontLeg.joints, rightFrontLeg.joints];
	var backLegJoints = [leftLegBack.joints, rightLegBack.joints]

	bodies.push(centralBody);

    return new Organism(genome, bodies, centralBody, joints, sensors);
}


