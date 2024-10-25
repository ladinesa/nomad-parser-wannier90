import os
import re
from typing import TYPE_CHECKING, Optional, List, Dict, Any

import numpy as np
from nomad.config import config

if TYPE_CHECKING:
    from nomad.datamodel import EntryArchive
    from structlog.stdlib import BoundLogger

from .parser import WOutParser
from .win_parser import WInParser
from nomad.parsing.file_parser.mapping_parser import TextParser, MetainfoParser
from nomad.parsing.file_parser.text_parser import DataTextParser
from nomad_parser_wannier90.schema_packages.package import Simulation, IN_ANNOTATION_KEY, OUT_ANNOTATION_KEY, BAND_ANNOTATION_KEY
from nomad_simulations.schema_packages.workflow import SinglePoint
from .utils import get_files, parse_dft_plus_tb_workflow

re_n = r'[\n\r]'

configuration = config.get_plugin_entry_point(
    'nomad_parser_wannier90.parsers:parser_entry_point'
)


class WBandTextParser(TextParser):

    def get_data(self, data: np.ndarray) -> np.ndarray:
        return np.transpose(data)[1:].transpose()



class WOutTextParser(TextParser):

    def get_lattice_vectors(self, vectors: List[Any]) -> np.ndarray:
        return np.vstack(vectors[-3:])

    def is_maximally_localized(self, niter: int, default=0) -> bool:
        return (niter or default) > 1



class WInTextParser(TextParser):

    def get_projections(self, source: List[Any]) -> List[Dict[str, Any]]:
        return [dict(atom=val) for val in source]

    def get_branch_label_indices(self, atom: List[Any], positions: List[np.ndarray], labels: List[str], lattice_vectors: List[np.ndarray]) -> Any:

        symbols, indices = [], []
        if not atom:
            return None

        elif isinstance(atom[0], int):
            indices = atom

        elif match := re.match(r'([cf])=(.+?),(.+?),(.+)', atom[0]):
            coord = match.groups()[0]
            position = np.array(match.groups()[1:4], float)
            if coord.lower() == 'f':
                position = np.dot(position, lattice_vectors)
            for n, pos in enumerate(positions):
                if np.allclose(position, pos, configuration.equal_cell_positions_tolerance):
                    indices.append(n)
                    symbols.append(labels[n])

        return dict(label=''.join(symbols), indices=indices)


class Wannier90Parser:
    def get_dft_files(self, mainfile:str) -> List[str]:
        for filename in ['vasprun.xml', 'OUTCAR']:
            files = get_files(
                pattern=f'*{filename}',
                filepath=mainfile,
                stripname=os.path.basename(mainfile),
                deep=False
            )
            if files:
                 return files
        return []

    def get_mainfile_keys(self, **kwargs):
        """
        Generates extra `child_archives` to create the DFT+TB workflow if the conditions in `workflow_dft_files` are met.
        """
        dft_files = self.get_dft_files(kwargs.get('filename', ''))
        return ['DFTPlusTB_workflow'] if dft_files else True

    def parse(self, mainfile: str, archive: 'EntryArchive', logger: 'BoundLogger', child_archives: Dict[str, 'EntryArchive'] = {}) -> None:
        # define mapping parser interface to OutParser
        wout_parser = WOutTextParser(text_parser=WOutParser())
        wout_parser.filepath = mainfile

        # construct metainfo parser
        data = Simulation()
        data_parser = MetainfoParser()
        data_parser.annotation_key = OUT_ANNOTATION_KEY
        data_parser.data_object = data

        wout_parser.convert(data_parser)
        archive.data = data

        # parse input file
        win_parser = WInTextParser(text_parser=WInParser())
        if data.model_system:
            win_files = get_files(
                pattern='*.win', filepath=mainfile, stripname=os.path.basename(mainfile)
            )
            if len(win_files) > 1:
                logger.warning(
                    'Multiple `*.win` files found. We will parse the first one.'
                )
            if win_files is not None:
                win_parser.filepath = win_files[0]
                # need data from out
                for key in ['structure', 'lattice_vectors']:
                    win_parser.data[key] = wout_parser.data.get(key)
                data_parser.annotation_key = IN_ANNOTATION_KEY
                data_parser.data_object = data
                win_parser.convert(data_parser)

        wband_parser = WBandTextParser(text_parser=DataTextParser())
        # parse band file
        band_files = get_files(
            pattern='*band.dat', filepath=mainfile, stripname=os.path.basename(mainfile)
        )
        for band_file in band_files:
            wband_parser.filepath = band_file
            wband_parser.data_object.parse('data')
            data_parser.annotation_key = BAND_ANNOTATION_KEY
            data_parser.data_object = data
            wband_parser.convert(data_parser)

        workflow = SinglePoint()
        workflow.normalize(archive=archive, logger=logger)
        archive.workflow2 = workflow

        # workflow
        if child_archives:
            from nomad.app.v1.routers.uploads import get_upload_with_read_access
            from nomad.datamodel import User

            upload_id = archive.metadata.upload_id
            upload = get_upload_with_read_access(
                upload_id=upload_id,
                user=User.get(user_id=archive.metadata.main_author.user_id),
            )
            dft_archive = None
            dft_files = self.get_dft_files(mainfile)
            dft_path = dft_files[-1].split('raw/')[-1]
            with upload.entries_metadata() as entries_metadata:
                for metadata in entries_metadata:
                    if metadata.mainfile == dft_path:
                        dft_archive = upload.get_entry(metadata.entry_id)._parser_results
                        break
            dft_plus_tb_archive = child_archives.get(
                'DFTPlusTB_workflow'
            )
            dft_plus_tb = parse_dft_plus_tb_workflow(
                dft_archive=dft_archive, tb_archive=archive
            )
            dft_plus_tb_archive.workflow2 = dft_plus_tb

        # debug
        self.wout_parser = wout_parser
        self.data_parser = data_parser
        self.win_parser = win_parser
        self.wband_parser = wband_parser
        # close parser contexts
        # wout_parser.close()
        # data_parser.close()
