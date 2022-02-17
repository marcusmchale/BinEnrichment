from typing import Optional, Set

from enum import Enum
from collections import defaultdict
import csv


class Expression(Enum):
	UP = 'up'
	DOWN = 'down'
	DIFF = 'diff'
	UNDETERMINED = 'undetermined'


class ExpressionMap:
	def __init__(self):
		self.undetermined = set()
		self.up = set()
		self.down = set()
		self.map = {
			'undetermined': self.undetermined,
			'up': self.up,
			'down': self.down,
			'diff': self.diff,
			'detected': self.detected
		}

	@property
	def detected(self):
		return set.union(self.undetermined, self.up, self.down)

	@property
	def diff(self):
		return set.union(self.up, self.down)


class Node:
	def __init__(
			self,
			bin_code,
			bin_name,
			bin_identifier,
			bin_description,
			parent
	):
		self.code: str = bin_code
		self.name: str = bin_name
		self.identifier: str = bin_identifier
		self.description: str = bin_description

		self.parent: Node = parent
		if parent:
			parent.add_child(self)
		self.children = set()

		self.genes = set()
		self.expression_map = ExpressionMap()

	def __hash__(self):
		return hash(self.code)

	def add_child(self, child: 'Node'):
		self.children.add(child)

	def add_gene(self, gene: str):
		self.genes.add(gene)

	def add_expression(self, gene: str, expression: Expression):
		self.expression_map.map[expression.value].add(gene)
		if self.parent:
			self.parent.add_expression(gene, expression)


class Tree:
	def __init__(self):
		self.root_code = "0"
		self.root = Node(self.root_code, "root", "root", "root", None)
		self.code_to_node = dict()
		self.gene_to_codes = defaultdict(set)
		self.unmapped_genes = set()

	def load_map(self, mapping_file_path):
		self.code_to_node[self.root_code] = self.root
		print('Start reading input file from: ' + mapping_file_path)
		with open(mapping_file_path) as tsv_file:
			reader = csv.reader(tsv_file, delimiter='\t', quoting=csv.QUOTE_NONE)
			next(reader)  # skip header, i.e. BINCODE	NAME	IDENTIFIER	DESCRIPTION	TYPE
			for row in reader:
				code = row[0].replace("'", "")
				name = row[1].replace("'", "")
				gene = row[2].replace("'", "")
				bin_description = row[3].replace("'", "")
				node: Optional[Node] = self.code_to_node.get(code)
				if not node:
					parent_code = code[0:code.rfind('.')] if '.' in code else self.root_code
					node = Node(
						code,
						name,
						gene,
						bin_description,
						self.code_to_node[parent_code]
					)
					self.code_to_node[code] = node
				if gene:
					node.add_gene(gene)
					self.gene_to_codes[gene].add(code)

	def add_degs(self, genes: Set[str], expression: Expression):
		for gene in genes:
			gene = gene.casefold()
			if gene not in self.gene_to_codes:
				self.unmapped_genes.add(gene)
			else:
				codes = self.gene_to_codes[gene]
				for code in codes:
					self.code_to_node[code].add_expression(gene, expression)

