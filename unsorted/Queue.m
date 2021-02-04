classdef Queue < handle
%QUEUE Time efficient Queue for fixed-length tuples of doubles
%
%	Version: 1.0
%	Date: 11/06/2020
%	Author: Thomas Thayer (tthayer@ucmerced.edu)
%
%   This class handles using a queue where every element in the queue has a fixed length and is of type 'double'.
	
	properties (Access = private)
		topIdx
		addBlock
		addIdx
		blockSize
		blocks
	end
	
	methods (Access = public)
		function obj = Queue(blockSize)
		%QUEUE Construct an instance of the Queue
		%   Inputs:
		%		blockSize: a vector of length 2 that defines how much memory to allocate for each block
		%			blockSize(1): number of queue elements in each block (default = 10000)
		%			blockSize(2): number of items in each element (default = 1)
		%	Outputs:
		%		obj: the newly constructed queue
			if ~exist('blockSize', 'var') || isempty(blockSize)
				blockSize = [10000, 1];
			elseif size(size(blockSize), 2) > 2
				error('blockSize must not contain more than two values');
			end
			obj.blockSize = blockSize;
			obj.topIdx = 1;
			obj.addBlock = 1;
			obj.addIdx = 1;
			obj.blocks = {zeros(blockSize)};
		end
		function add(obj, element)
		%ADD Add an element or set of elements to the queue
		%   Inputs:
		%		element: a matrix containing the elements to add to the queue, where rows contain individual elements and columns contain the items for each element
			if ~exist('element', 'var') || isempty(element)
				warning('Nothing to add to queue.');
			end
			if (size(element, 2) ~= obj.blockSize(2))
				error('New elements must be of appropriate size.');
			end
			while ~isempty(element)
				spaceLeft = obj.blockSize(1) - obj.addIdx + 1;
				if size(element, 1) <= spaceLeft
					tempSet = element;
					element = [];
				else
					tempSet = element(1:spaceLeft, :);
					element(1:spaceLeft, :) = [];
				end
				obj.blocks{obj.addBlock}(obj.addIdx:(obj.addIdx+size(tempSet, 1)-1), :) = tempSet;
				obj.addIdx = obj.addIdx + size(tempSet, 1);
				if obj.addIdx > obj.blockSize(1)
					obj.addBlock = obj.addBlock + 1;
					obj.addIdx = 1;
					obj.blocks{obj.addBlock} = zeros(obj.blockSize);
				end
			end
		end
		function topElement = top(obj)
		%TOP Get the next item in the queue
		%	Outputs:
		%		topElement: the item in position 1 of the queue
			if obj.isEmpty()
				warning('Queue is empty.');
				topElement = [];
			else
				topElement = obj.blocks{1}(obj.topIdx, :);
			end
		end
		function topElement = pop(obj)
		%POP Get and remove the next item in the queue
		%	Outputs:
		%		topElement: the item removed from the queue, formerly at position 1
			topElement = top(obj);
			if ~obj.isEmpty()
				if obj.topIdx == obj.blockSize(1)
					obj.topIdx = 1;
					if length(obj.blocks) > 1
						obj.blocks(1) = [];
						obj.addBlock = obj.addBlock-1;
					end
				else
					obj.topIdx = obj.topIdx + 1;
				end
			else
				error('Cannot pop queue. Queue is empty.');
			end
		end
		function blockSize = getBlockSize(obj)
		%GETBLOCKSIZE Return the size of each block
		%	Outputs:
		%		blockSize: a vector of length 2 that defines how much memory to allocate for each block
		%			blockSize(1): number of queue elements in each block (default = 10000)
		%			blockSize(2): number of items in each element (default = 1)
			blockSize = obj.getBlockSize;
		end
		function maybe = isEmpty(obj)
		%ISEMPTY Is the queue empty
		%	Outputs:
		%		maybe: a boolean indicating if the queue is empty
			if (obj.addBlock == 1) && (obj.addIdx == obj.topIdx)
				maybe = true;
			else
				maybe = false;
			end
		end
		function queueList = getList(obj)
		%GETLIST Get a list of items in the queue
		%	Outputs:
		%		queueList: a list containing the elements in the queue, where rows contain individual elements and columns contain the items for each element
			queueList = zeros(obj.getSize());
			counter = 1;
			for i=1:obj.addBlock
				if i==1
					this_size = size(obj.blocks{i}(obj.topIdx:end, :), 1);
					queueList(counter:(counter+this_size-1), :) = obj.blocks{i}(obj.topIdx:end, :);
					counter = counter + this_size;
				elseif i==obj.addBlock
					this_size = size(obj.blocks{i}(1:(obj.addIdx-1), :), 1);
					queueList(counter:(counter+this_size-1), :) = obj.blocks{i}(1:(obj.addIdx-1), :);
					counter = counter + this_size;
				else
					this_size = size(obj.blocks{i}, 1);
					queueList(counter:(counter+this_size-1), :) = obj.blocks{i};
					counter = counter + this_size;
				end
			end
		end
		function queueSize = getSize(obj)
		%GETSIZE Get the number of elements in the queue
		%	Outputs:
		%		queueSize: an integer indicating how many elements are in the queue
			queueSize = [obj.blockSize(1)*(obj.addBlock - 1) + obj.addIdx - obj.topIdx, obj.blockSize(2)];
		end
	end
end

