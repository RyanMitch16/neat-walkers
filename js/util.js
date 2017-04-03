
var HashMap = function() {
	var items = [undefined];
	var itemCount = 0;

	/** Retrieve the count of key-value pairs stored in the hash map. */
	this.getSize = function() {
		return itemCount;
	}

	this.isEmpty = function() {
		return itemCount == 0;
	}

	this.getEntries = function() {
		var output = [];
		for (var i = 0; i < items.length; i++) {
			if (items[i] != undefined) {
				output.push(items[i]);
			}
		}
		return output;
	}

	this.get = function(key) {
		return items[find(key)]["value"];
	}

	this.clear = function() {
		items = [undefined];
		itemCount = 0;
	}

	/** */
	this.put = function(key, value) {
		if (this.contains(key)) {
			throw "The key " + key + " has previously been added to the hash map";
		}

		// Grow the size of the array and copy the key-value pairs back into it
		if (itemCount == items.length) {
			var arr = [];
			for (var i = 0; i < items.length * 2; i++) {
				arr[i] = undefined;
			}
			for (var i = 0; i < items.length; i++) {
				insert(arr, items[i]["key"], items[i]["value"]);
			}
			items = arr;
		}
		insert(items, key, value);
		itemCount += 1;
	}

	/** */
	this.remove = function(key) {
		if (!this.contains(key)) {
			throw "The key " + key + " is not present in the hash map"
		}

		// Shrink the size of the array and copy the key-value pairs back into it
		if (itemCount * 2 == items.length) {
			var arr = [];
			for (var i = 0; i < items.length / 2; i++) {
				arr[i] = undefined;
			}

			for (var i = 0; i < items.length; i++) {
				if (items[i] != undefined) {
					insert(arr, items[i]["key"], items[i]["value"]);
				}
			}
			items = arr;
		}

		items[find(key)] = undefined;
		itemCount -= 1;
	}

	/** Checks if a key-value pair with the specified key exists in the hash map already. */
	this.contains = function(key) {
		return find(key) != -1;
	}

	/** */
	var find = function(key) {
		var hashCode = computeHashCode(key);
		var i = 0;
		for (var i = 0; i < items.length; i++) {
			pos = (hashCode + i) % items.length;
			if (items[pos] != undefined && checkEquality(items[pos]["key"], key)) {
				return pos;
			}
		}
		return -1;
	}

	/** */
	var insert = function(arr, key, value) {
		// Find a free location for the key using linear probing
		var hashCode = computeHashCode(key);
		var i = 0;
		do {
			pos = (hashCode + i) % arr.length;
			i += 1;
			//console.log(pos+" "+arr.length);
		} while (arr[pos] != undefined);

		// Insert a map of the key-value pair into the array
		arr[pos] = {"key" : key, "value" : value};
	}

	/** Performs a depp equality check of the specified values. */
	var checkEquality = function(obj1, obj2) {
		if (obj1 === obj2 || obj1 == obj2) {
			return true;
		} else if (obj1 == undefined || obj1 == null || obj2 == undefined || obj2 == null) {
			return false;
		} else if (obj1.equals != undefined && obj2.equals != undefined) {
			return obj1.equals(obj2) && obj2.equals(obj1);
		}
		return false;
	} 

	/** Computes the integer has code of a number, string, or an object. */
	var computeHashCode = function(key, d) {
		if (d == 50) return 0;
		if (d == undefined) d = 0;

		// Compute the hash value of a a number
		if (typeof key === 'number') {
			return key;
		}
		// Compute the hash value of a string
		if (typeof key === 'string') {
			var hash = 7;
			for (var i = 0; i < key.length; i++) {
			    hash = hash * 31 + key.charCodeAt(i);
			}
			return hash;
		}
		if (typeof key === 'object' && key != null && key != undefined) {
			// Check if the object already has a hash code function
			if (key.hasOwnProperty("hashCode")) {
				hashCode = key.hashCode();
			}
			// Compute the hash value of an object
			var hash = 7;
			var keys = Object.keys(key).sort();
			//console.log(key)
			for (var i = 0; i < keys.length; i++) {
			    hash = hash * 31 & computeHashCode(key[keys[i], d + 1]);
			}
			return hash;
		}
		return 0;
	}
}