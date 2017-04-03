

/** Randomly shuffle the order of the elements in the provided array. */
var randomizeArrayOrder = function(a) {
	for (var i = 0; i < a.length * 7; i++) {
		var j = Math.floor(Math.random() * a.length);
		var k = Math.floor(Math.random() * a.length);
		var temp = a[j];
		a[j] = a[k];
		a[k] = temp;
	}
	return a;
}

/** */
var NeuralNetwork = function(inputCount, outputCount) {

	/** Initializes a new neuron. */
	class Neuron {
		constructor(id){
			this.id = id;
			this.value = 0;
			this.layer = 0;
			this.inputNeurons = new HashMap();
			this.outputNeurons = new HashMap();
		}
	}

	/** Updates the layer values of neurons. */
	var updateLayersOfFrom = function(n) {
		var outputNeurons = n.outputNeurons.getEntries();
		for (var i = 0; i < outputNeurons.length; i++) {
			if (n.layer >= outputNeurons[i]["value"].layer) {
				outputNeurons[i]["value"].layer = n.layer + 1;
				updateLayersOfFrom(outputNeurons[i]["value"]);
			}
		}
	}

	// The list of neurons in the network indexed by id
	var neurons = new HashMap();

	// The lists of lists of neurons representing the order in which neuron values must be calculate to ensure neuron
	// dependencies are computed in order 
	var cachedLayers = undefined;

	// Initialize the input and output neurons for the network
	for (var i = 0; i < inputCount + outputCount; i++) {
		neurons.put(i, new Neuron(i));
	}

	/** Returns the neuron with the specified id if it exists in the neural network. */
	var getNeuron = function(id) {
		if (neurons.contains(id)) {
			return neurons.get(id);
		} else {
			throw "Neuron with id " + id + " does not exists";
		}
	}

	/** */
    var checkCycle = function(currentNeuron, target) {
      var outputNeurons = currentNeuron.outputNeurons.getEntries();
      for (var i = 0; i < outputNeurons.length; i++) {
         if (target == outputNeurons[i]["key"]) {
            return true;
         } else {
            if (checkCycle(outputNeurons[i]["value"], target)) {
            	return true;
            }
         }
      }
      return false;
   }

	this.containsNeuron = function(id) {
		return neurons.contains(id);
	}

	/** Creates a new neuron that will be used in the neural network. */
	this.addNeuron = function(id) {
		neurons.put(id, new Neuron(id));
	}

	/** Adds a connection between the nodes of specified ids with the specified edge weight. */
	this.addConnection = function(from, to, weight) {
		// Check if the to node is an input node
		if (to < inputCount) {
			throw "Input nodes cannot have incoming connections " + to;
		}
		// Checks if the connection already has been made
		if (this.hasConnection(from, to) || this.hasConnection(to, from)) {
			throw "Connection from " + from + " to " + to + " already exists";
		}

		// Add the connections to the hashmaps of both neurons
		getNeuron(from).outputNeurons.put(to, getNeuron(to));
		getNeuron(to).inputNeurons.put(from, {"neuron" : getNeuron(from), "weight" : weight, "enabled" : true});

		// Clear the cached layer or neurons because the dependency order may have changed
		cachedLayers = undefined;
		updateLayersOfFrom(getNeuron(from));
	}

   /** */
   this.isConnectionLegal = function(from, to) {
      return !this.containsNeuron(from) 
      	|| !this.containsNeuron(to) 
      	|| ((from < inputCount || from >= inputCount + outputCount)
         	&& (to >= inputCount)
         	&& (from != to)
         	&& !this.hasConnection(from, to)
         	&& !checkCycle(getNeuron(to), from));
   }

	/** Checks if a connection between the neurons of the specified ids exists or not. */
	this.hasConnection = function(from, to) {
		return getNeuron(to).inputNeurons.contains(from);
	}

	/** Retrieves the weight of the connection between the neurons of the specified ids. */
	this.setConnectionWeight = function(from, to, weight) {
		if (this.hasConnection(from, to)) {
			getNeuron(to).inputNeurons.get(from)["weight"] = weight;
		} else {
			throw "Connection from " + from + " to " + to + " does not already exists";
		}
	}

	/** */
	this.removeConnection = function(from, to) {
		if (this.hasConnection(from, to)) {
			getNeuron(to).inputNeurons.remove(from);
			getNeuron(from).outputNeurons.remove(to);
		} else {
			throw "Connection from " + from + " to " + to + " does not already exists";
		}
	}

	/** Retrieves the identifiers of the input neurons. */
	this.getInputIds = function() {
		var ins = [];
		for (var i = 0; i < inputCount; i++) {
			ins.push(i);
		}
		return ins;
	}

	/** Retrieves the identifiers of the output neurons. */
	this.getOutputIds = function() {
		var outs = [];
		for (var i = 0; i < outputCount; i++) {
			outs.push(i + inputCount);
		}
		return outs;
	}

	this.getNeuronIds = function() {
		var neuronEntries = neurons.getEntries();
		var outs = [];
		for (var i = 0; i < neuronEntries.length; i++) {
			outs.push(neuronEntries[i]["key"]);
		}
		return outs;
	}

	/** 
	 * Proogates the input values through the network.
	 * @param {number[]} inputValues - the in-order values of the input nodes of the network
	 * @returns {number[]} the values of the output nodes after propogation
	 */
	this.propagate = function(inputValues) {
		var neuronEntries = neurons.getEntries();

		// Set the input values of the network
		if (inputValues.length != inputCount) {
			throw "Incorrect number of input values provided (" + inputValues.length + " should be " + inputCount + ")";
		}

		for (var i = 0; i < inputCount; i++) {
			getNeuron(i).value = inputValues[i];
		}

		// Sets up the mapping of layer number to the list of nodes contained in that layer
		if (cachedLayers == undefined) {
			cachedLayers = [];
			for (var i = 0; i < neuronEntries.length; i++) {
				var n = neuronEntries[i]["value"];
				if (cachedLayers[n.layer] == undefined) {
					cachedLayers[n.layer] = [];
				}
				cachedLayers[n.layer].push(n);
			}
		}

		// For each layer, compute the value of each node from the values of its imput node
		for (var i=0; i<cachedLayers.length; i++) {
			for (var j=0; j<cachedLayers[i].length; j++) {
				neuron = cachedLayers[i][j];
				if (neuron.id < inputCount) {
					continue;
				}

				// Calculate the total weight of all the inputs
				var x = 0;
				var inputNeurons = neuron.inputNeurons.getEntries();
				for (var k = 0; k < inputNeurons.length; k++) {
					var input = inputNeurons[k]["value"];
					if (input["enabled"]) {
						x += input.neuron.value * input["weight"];
					}
				}

				//x  = Math.random() * 2 - 1;

				// Use the sigmoid function as the activation function
				neuron.value = 1 / (1 + Math.pow(Math.E, -x));
			}
		}

		// Return the values from the output layer
		var output = [];
		for (var i = inputCount; i < inputCount + outputCount; i++) {
			output.push(getNeuron(i).value);
		}
		return output;
	}

	/** Creates a deep copy of the neural network. */
	/*this.clone = function() {
		var nn = new NeuralNetwork(inputCount, outputCount);

		// Add the hidden nodes to the neural network
		var neuronEntries = neurons.getEntries();
		for (var i = 0; i < neuronEntries.length; i++) {
			if (neuronEntries[i]["key"] >= inputCount + outputCount) {
				nn.addNeuron(neuronEntries[i]["key"]);
			}
		}

		// Add the input connections for all the non-input neurons
		for (var i = 0; i < neuronEntries.length; i++) {
			var inputNeurons = neuronEntries[i]["value"].inputNeurons.getEntries();
			for (var k = 0; k < inputNeurons.length; k++) {
				nn.addConnection(inputNeurons[k]["key"], neuronEntries[i]["key"], inputNeurons[k]["value"]["weight"]);
			}
		}
		return nn;
	}*/
}

/** */
var EdgeGene = function(from, to, weight, enabled) {
	this.from = from;
	this.to = to;
	this.weight = weight;
	this.enabled = enabled;
}

class InnovationNumber {

	/** */
	constructor(initial) {
		this._number = initial;
	}

	/** */
	next() {
		return this._number++;
	}
}

/** */
class Genome {

	/** Creates a new Genome with the specified number of inputs and outputs the genome neural network requires. */
	constructor(inputCount, outputCount, innovationNumber) {
		this._innovationNumber = innovationNumber;
		this._inputCount = inputCount;
		this._outputCount = outputCount;
		this._genome = new HashMap();
		this._values = [];
		this._network = new NeuralNetwork(inputCount, outputCount);
	}

	/** 
	 * Retrieves the list identifiers of the input neurons of underlying neural network of the genome.
	 * @returns {number[]}
	 */
	getInputIds() {
		return this._network.getInputIds();
	}

	/**
	 * Retrieves the list identifiers of the output neurons of the underlying neural network of the genome. 
	 * @returns {number[]}
	 */
	getOutputIds() {
		return this._network.getOutputIds();
	}

	/** 
	 * Adds a new gene to the Genome representing a connection in the neural network between two neurons. If adding the 
	 * connection would cause a cycle in the neural network, then the connection is disabled when added.
	 * @param from {number} the identifier of the neuron the connection is coming from
	 * @param to {number} the identifier of the neuron the connection is going to
	 * @param weight {number} the weight of the ceonnection edge
	 * @param enabled {boolean} whether the connection should be enabled or not
	 */
	addEdgeGene(from, to, weight, enabled) {
		this._addEdgeGene(this._innovationNumber.next(), from, to, weight, enabled);
	}

	/** 
	 * The private implementation of adding a connection gene that allows specification of the innovation number to 
	 * give the gene.
	 */
	_addEdgeGene(innovationNumber, from, to, weight, enabled) {
		// Add the neurons to the networ if they do not already exist
		if (!this._network.containsNeuron(from)) {
			this._network.addNeuron(from);
		}
		if (!this._network.containsNeuron(to)) {
			this._network.addNeuron(to);
		}

		// Add the connection only if the gene is enabled
		if (enabled) {
			if (this._network.isConnectionLegal(from, to)) {
				this._network.addConnection(from, to, weight);
			} else {
				enabled = false;
			}
		}

		// Add the gene to the gene hashmap
		this._genome.put(innovationNumber, new EdgeGene(from, to, weight, enabled));
	}

	/** 
	 * Adds a new gene to the Genome representing a variable that does not represnt a connection and that must be fit
	 * to a value.
	 * @param max {number} the maximum value of the variable
	 * @param min {number} the minimum value of the variable
	 */
	addVariableGene(min, max) {
		this._values.push({
			"value" : (max - min) / 2 + min,
			"min" : min,
			"max" : max
		})
	}

	/** Retrieves the value of the nth variable gene added to the Genome. */
	getVariableGene(n) {
		return this._values[n].value;
	}

	/** Retrives a deep copy of the genome. */
	clone() {
		var newGenome = new Genome(this._inputCount, this._outputCount, this._innovationNumber);

		// Add the values of the edge genes to the copy
		var genomeEntries = this._genome.getEntries();
		for (var i = 0; i < genomeEntries.length; i++) {
			var gene = genomeEntries[i].value;
			newGenome._addEdgeGene(genomeEntries[i].key, gene.from, gene.to, gene.weight, gene.enabled);
		}

		// Add the values of the value genes to the copy
		for (var i = 0; i < this._values.length; i++) {
			newGenome.addVariableGene(this._values[i].value, this._values[i].min, this._values[i].max);
		}
		return newGenome;
	}

	/** 
	 * Propogates the values of the input neurons through the neural network.
	 * @param inputValues {number[]} the values of the input neurons
	 * @returns {number[]} the values of the output neurons after the propogation of the input values
	 */
	propagate(inputValues) {
		return this._network.propagate(inputValues);
	}

	/** Mutates the weights of the edges and variables of the Genomes. */
	mutate() {
		this._mutateVariables();
		this._mutateWeights();
		if (Math.random() < 0.05) {
			this._mutateAddConnection();
		}
		if (Math.random() < 0.03) {
			this._mutateAddNeuron();
		}
	}

	/** Adds a new connection between two currently unconnected neurons in the network. */
	_mutateAddConnection() {
		var neurons = this._network.getNeuronIds();
		var genomeEntries = this._genome.getEntries();
		var suffledNeurons1 = randomizeArrayOrder(neurons.slice());
		var suffledNeurons2 = randomizeArrayOrder(neurons.slice());
		for (var i = 0; i < suffledNeurons1.length; i++) {
			for (var j = 0; j < suffledNeurons2.length; j++) {
				if (this._network.isConnectionLegal(suffledNeurons1[i], suffledNeurons2[j])) {
					// Check if the gene is not present and currently disabled
					var foundDisabled = false;
					for (var k = 0; k < genomeEntries.length; k++) {
						var gene = genomeEntries[k].value;
						if (gene.from == suffledNeurons1[i] && gene.to == suffledNeurons2[j]) {
							foundDisabled = true;
							break;
						}
					}

					// Add the gene with a random weight if the connection does not currently exist
					if (!foundDisabled) {
						var weight = Math.random() * 2 - 1;
						this.addEdgeGene(suffledNeurons1[i], suffledNeurons2[j], weight, true);
						return;
					}
				}
			}
		}
	}

	/** Splices a current connection by inserting a new neuron between two currently connected nodes. */
	_mutateAddNeuron() {
		var genomeEntries = this._genome.getEntries();
		if (genomeEntries.length > 0) {
			// Find an enabled connection to mutate
			var gene;
			do {
				gene = genomeEntries[Math.floor(Math.random() * genomeEntries.length)].value;
			} while (!gene.enabled);

			// Disable the original connection
			gene.enabled = false;

			// Add the connection from the from neuron to the new neuron
			var id = this._innovationNumber.next() + this._inputCount + this._outputCount;
			this.addEdgeGene(gene.from, id, 1, true);
			this.addEdgeGene(id, gene.to, gene.weight, true);
		}
	}

	/** Mutates the values of the Genome variables. */
	_mutateVariables() {
		for (var i = 0; i < this._values.length; i++) {
			if (Math.random() < 0.9) {
				var range = this._values[i].max - this._values[i].min;
				if (Math.random() < 0.9) {
					// Clamp the perturbed value to the specified range
					var newValue = this._values[i].value + (Math.random() * (range / 50)) - (range / 25);
					this._values[i].value = Math.min(Math.max(newValue, this._values[i].min), this._values[i].max);
				} else {
					// Randomly assign a new value within the specified range
					this._values[i].value = (Math.random() * range) + this._values[i].min;
				}
			}
		}
	}

	/** Randomly mutate the enabled weights of the genome. */
	_mutateWeights() {
		var genomeEntries = this._genome.getEntries();
		for (var i = 0; i < genomeEntries.length; i++) {
			var gene = genomeEntries[i].value;
			if (gene.enabled) {
				if (Math.random() < 0.9) {
					if (Math.random() < 0.9) {
						// Randomly perturbe the gene weight
						gene.weight = gene.weight + Math.random() * 0.2 - 0.1;
					} else {
						// Randomly assign a new gene weight
						gene.weight = Math.random() * 4 - 2;
					}

					// Set the weight of the connection in the network
					this._network.setConnectionWeight(gene.from, gene.to, gene.weight);
				}
			}
		}
	}

	/** Creates a offsping Genome from the mating of this dominant genome with the specified recessive genome. */
	crossOver(recessive) {
		var dominantGenomeEntries = this._genome.getEntries().sort(function(a, b){return a.key - b.key});
		var recessiveGenomeEntries = recessive._genome.getEntries().sort(function(a, b){return a.key - b.key});
		var newGenome = new Genome(this._inputCount, this._outputCount, this._innovationNumber);

		// Randomly choose which genome variable to use for each variable
		for (var i = 0; i < this._values.length; i++) {
			var variable = (Math.random() >= 0.5) ? this._values[i] : recessive._values[i];
			newGenome.addVariableGene(variable.value, variable.min, variable.max);
		}

		var d = 0, r = 0;
		while (d < dominantGenomeEntries.length || r < recessiveGenomeEntries.length) {
			// There are no more excess genes of the dominant genome to add
			if (d >= dominantGenomeEntries.length) {
				return newGenome;
			} 
			// Add excess genes from the dominant genome
			else if (r >= recessiveGenomeEntries.length) {
				var gene = dominantGenomeEntries[d].value;
				var enabled = gene.enabled ? true : ((Math.random() < 0.25) ? true : false);
				newGenome._addEdgeGene(dominantGenomeEntries[d].key, gene.from, gene.to, gene.weight, enabled);
				d += 1;
			} else {
				// Randomly assign non-disjoint and non-excess genes
				if (dominantGenomeEntries[d].key == recessiveGenomeEntries[r].key) {
					var gene = (Math.random() >= 0.5) ? dominantGenomeEntries[d].value : recessiveGenomeEntries[r].value;
					var enabled = gene.enabled ? true : ((Math.random() < 0.25) ? true : false);
					newGenome._addEdgeGene(dominantGenomeEntries[d].key, gene.from, gene.to, gene.weight, enabled);
					d += 1;
					r += 1;
				// Add disjoint genes from the dominant genome
				} else if (dominantGenomeEntries[d].key < recessiveGenomeEntries[r].key) {
					var gene = dominantGenomeEntries[d].value;
					var enabled = gene.enabled ? true : ((Math.random() < 0.25) ? true : false);
					newGenome._addEdgeGene(dominantGenomeEntries[d].key, gene.from, gene.to, gene.weight, enabled);
					d += 1;
				// Do not add disjoint genes from the recessive genome
				} else {
					r+= 1;
				}
			}
		}

		return newGenome;
	}

	/** Checks if this genome and the specified other genome are close enough to be considered the same species. */
	isSameSpecies(other) {
		var dominantGenomeEntries = this._genome.getEntries().sort(function(a, b){return a.key - b.key});
		var recessiveGenomeEntries = other._genome.getEntries().sort(function(a, b){return a.key - b.key});

		var matched = 0, averageWeightDifference = 0, disjoint = 0, excess = 0, d = 0, r = 0;
		while (d < dominantGenomeEntries.length || r < recessiveGenomeEntries.length) {
			// Add excess genes from the dominant genome
			if (d >= dominantGenomeEntries.length) {
				excess += 1;
				r+= 1;
			} 
			// Add excess genes from the recessive genome
			else if (r >= recessiveGenomeEntries.length) {
				excess += 1;
				d += 1;
			} else {
				// Add dominant genes over recessive genes
				if (dominantGenomeEntries[d].key == recessiveGenomeEntries[r].key) {
					matched += 1;
					averageWeightDifference += Math.abs(dominantGenomeEntries[d].value.weight - 
						recessiveGenomeEntries[r].value.weight);
					d += 1;
					r += 1;
				// Add disjoint genes from the dominant genome
				} else if (dominantGenomeEntries[d].key < recessiveGenomeEntries[r].key) {
					disjoint += 1;
					d += 1;
				// Add disjoint genes from the recessive genome
				} else {
					disjoint += 1; 
					r+= 1;
				}
			}
		}

		var c1 = 1, c2 = 1, c3 = 0.4, N = 5, threshold = 3;
		if (matched > 0 ) {
			averageWeightDifference = averageWeightDifference / matched;
		}
		dominantGenomeEntries = null;
		recessiveGenomeEntries = null;
		return ((((c1 * excess) + (c2 * disjoint)) /  N) + (c3 * averageWeightDifference))
			< threshold;
	}
}

var GenomePool = function(size, resetFunction) {

	var Species = function(representativeGenome) {
		this.genomes = new HashMap();
		this.maxFitness = 0;
		this.stagnemntGenerations = 0;

		this.genomes.put(representativeGenome, representativeGenome);
		this.representativeGenome = representativeGenome;
	}

	var species = [];
	var fitnesses = new HashMap();

	var resetFunction = resetFunction;

	resetFunction();

	this.getTotalSpecies = function() {
      return species.length;
   }

	this.assignFitness = function(genome, fitness) {
		if (fitnesses.contains(genome)) {
			console.log(genome === fitnesses.get(genome));

			throw "Fitness for the genome " + genome + " already assigned.";
		} else {
			fitnesses.put(genome, fitness);
		}
	}

	this.startTesting = function() {
		fitnesses.clear();
	}

	this.endTesting = function() {
		/*if (fitnesses.size() != genomes.size()) {
			throw "Fitness for " + (genomes.size() - fitnesses.size()) +" genomes are missing.";
		}*/

		var genomeFitnessEntries = fitnesses.getEntries();

		for (var i = 0; i < genomeFitnessEntries.length; i++) {
			var genome = genomeFitnessEntries[i]["key"];

			// Determine the species the genome belongs to 
			for (var j = 0; j < species.length; j++) {
				if (species[j].genomes.contains(genome)) {
					break;
				} 
				if (genome.isSameSpecies(species[j].representativeGenome)) {
					species[j].genomes.put(genome, genome);
					break;
				}
			}
			// Make the genome the representative for its species
			if (j == species.length) {
				species.push(new Species(genome));
			}

			// 
			if (species[j].maxFitness < genomeFitnessEntries[i]["value"]) {
				species[j].maxFitness = genomeFitnessEntries[i]["value"];
				species[j].stagnemntGenerations = -1;
			}
		}

		// Increment the count of generations without improvement in maximum fitness
		for (var i = 0; i < species.length; i++) {
			species[i].stagnemntGenerations += 1;

			// Kill of species that are not improving
			if (species[i].stagnemntGenerations == 25) {
				species[i] = null;
			}
		}

		// Filter out the killed off species
		species = species.filter(function(a) { return a != null; });

		// Calculate adjusted fitnesses based on closeness to other organisms
		var sumOfAverageFitnesses = 0;
		var averageFitnesses = [];
		for (var i = 0; i < species.length; i++) {
			var genomeEntriess = species[i].genomes.getEntries();
			var sumFitness = 0;
			for (var j = 0; j < genomeEntriess.length; j++) {
				sumFitness = fitnesses.get(genomeEntriess[j]["value"]);
			}
			var avg = sumFitness / genomeEntriess.length;
			averageFitnesses.push(avg);
			sumOfAverageFitnesses += avg;
			genomeEntriess = null;
		}

		var bestPercentage = 0.50;
		var keepBestPercentage = 0.05;
		var nextGenerationGenomes = [];
		var bestGenomes = [];

		for (var i = 0; i < species.length; i++) {
			// Sort the genomes in descending order of fitness
			var sortedGenomes = species[i].genomes.getEntries().sort(function(a, b){
				return fitnesses.get(b["value"]) - fitnesses.get(a["value"]);
			});

			// Remove the less fit organisms from each species
			var keep = Math.min(Math.floor((size * bestPercentage) * (averageFitnesses[i] / sumOfAverageFitnesses)), 
				sortedGenomes.length)
			for (var j = 0; j < sortedGenomes.length; j++) {
				if (j < keep) {
					bestGenomes.push(sortedGenomes[j]["value"])
				} else {
					species[i].genomes.remove(sortedGenomes[j]["value"]);
				}
			}

			//console.log(species[i].genomes.getEntries().length);

			if (species[i].genomes.getSize() == 0) {
	            species[i] = null;
	         }

	         // Add the champion of each large species directly to the next generation
	         else if (species[i].genomes.getSize() >= Math.floor(keepBestPercentage * size)) {
	            nextGenerationGenomes.push(sortedGenomes[0]["value"]);
	         }
			sortedGenomes = null;
		}

		// Filter out the killed off species
		species = species.filter(function(a) { return a != null });

		if (species.length == 0) {
			console.log("RESET");
			fitnesses.clear();
			resetFunction();
			species = [];
			return null;
		}


		var onlyMutationPercentage = 0.25;


		console.log("Start")

		// Cross over and mutate organisms from 
		while (nextGenerationGenomes.length + Math.floor(onlyMutationPercentage * size) < size) {
			var genome1, genome2;
			if (Math.random() < 0.02) {
				// Provide a chance for inter species mating
				genome1 = bestGenomes[Math.floor((Math.random() * bestGenomes.length))];
				genome2 = bestGenomes[Math.floor((Math.random() * bestGenomes.length))];
			} else {
				// Crossover organisms of the same species
				var speciesGenomes = species[Math.floor(Math.random() * species.length)].genomes.getEntries();
				genome1 = speciesGenomes[Math.floor((Math.random() * speciesGenomes.length))]["value"];
				genome2 = speciesGenomes[Math.floor((Math.random() * speciesGenomes.length))]["value"];
				speciesGenomes = null;
			}

			console.log("s1")
			var resultingGenome;
			if (fitnesses.get(genome1) > fitnesses.get(genome2)) {
				resultingGenome = genome1.crossOver(genome2);
			} else if (fitnesses.get(genome1) < fitnesses.get(genome2)) {
				resultingGenome = genome2.crossOver(genome1);
			} else {
				resultingGenome = (Math.random() > 0.5) ? genome1.crossOver(genome2) : genome2.crossOver(genome1);
			}
			console.log("e1")


			resultingGenome.mutate();		
			nextGenerationGenomes.push(resultingGenome);
		}

		console.log("s2")

		while (nextGenerationGenomes.length < size) {
			var genome = bestGenomes[Math.floor((Math.random() * bestGenomes.length))].clone();
			genome.mutate();
			nextGenerationGenomes.push(genome);
		}

		console.log("e2")

		console.log("end")


		// 
		fitnesses.clear();

		for (var i = 0; i < species.length; i++) {
			species[i].genomes.clear();
		}


		return nextGenerationGenomes;

	}
}
