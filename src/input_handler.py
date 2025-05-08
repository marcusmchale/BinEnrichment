from pathlib import Path

class InputHandler:

    @staticmethod
    def read_file(file_path):
        with open(file_path, 'r') as f:
            return [l.rstrip() for l in f]

    def __init__(self, args):
        self.mapping = Path(args['mapping'])
        if not self.mapping.exists():
            raise ValueError("Mapping file does not exist")

        self.output = Path(args['output'])
        if self.output.exists():
            print('Output file exists and will be overwritten')

        self.background = args['background_list'] or []
        if args['background_file']:
            bf = Path(args['background_file'])
            if not bf.exists():
                raise ValueError("Background file not found")
            self.background += self.read_file(bf)

        if args['target_list'] or args['target_file']:
            self.target = args['target_list'] or []
            if args['target_file']:
                tf = Path(args['target_file'])
                if not tf.exists():
                    raise ValueError("Target file not found")
                self.target += self.read_file(tf)
        else:
            self.target = None

        if args['up_list'] or args['up_file']:
            self.up = args['up_list'] or []
            if args['up_file']:
                uf = Path(args['up_file'])
                if not uf.exists():
                    raise ValueError("Up file not found")
                self.up += self.read_file(uf)
        else:
            self.up = None

        if args['down_list'] or args['down_file']:
            self.down = args['down_list'] or []
            if args['down_file']:
                df = Path(args['down_file'])
                if not df.exists():
                    raise ValueError("Down file not found")
                self.down += self.read_file(df)
        else:
            self.down = None

        self.alpha = args['alpha']


