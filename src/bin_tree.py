import csv


class BinTree:
	def __init__(self):
		self.root_node = None
		self.bin_code_to_node = dict()
		self.target_to_bin_code = dict()
		self.root_bin_code = "0"

	def load_mapping_file(self, mapping_file_path):
		self.root_node = self.BinNode(self.root_bin_code, "root", "root", "root")
		self.bin_code_to_node[self.root_bin_code] = self.root_node
		print('Start reading input file from: ' + mapping_file_path)
		with open(mapping_file_path) as tsv_file:
			reader = csv.reader(tsv_file, delimiter='\t', quoting=csv.QUOTE_NONE)
			next(reader)  # skip header, i.e. BINCODE	NAME	IDENTIFIER	DESCRIPTION	TYPE
			for row in reader:
				bin_code = row[0].replace("'", "")
				bin_name = row[1].replace("'", "")
				target_identifier = row[2].replace("'", "")
				bin_description = row[3].replace("'", "")
				node = self.bin_code_to_node.get(bin_code)
				if not node:
					parent_code = bin_code[0:bin_code.rfind('.')] if '.' in bin_code else self.root_bin_code
					node = self.BinNode(
						bin_code,
						bin_name,
						target_identifier,
						bin_description,
						self.bin_code_to_node[parent_code]
					)
					self.bin_code_to_node[bin_code] = node
				if target_identifier:
					node.add_target(target_identifier)
					if target_identifier in self.target_to_bin_code:
						self.target_to_bin_code[target_identifier].append(node)
					else:
						self.target_to_bin_code[target_identifier] = [node]

	def get_nodes_for_target(self, target):
		return self.target_to_bin_code[target.lower()]

	class BinNode:
		def __init__(
				self,
				bin_code,
				bin_name,
				bin_identifier,
				bin_description,
				parent=None
		):
			self.code = bin_code
			self.name = bin_name
			self.identifier = bin_identifier
			self.description = bin_description
			self.children = []
			self.target_list = []
			self.parent = parent
			if parent:
				parent.add_child(self)

		def add_child(self, child):
			self.children.append(child)

		def add_target(self, target_name):
			self.target_list.append(target_name)
