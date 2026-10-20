# -*- coding: utf-8 -*-
import pandas as pd
import pytest

from automlsa2.phylogeny import check_outgroup

LABELS = ['Kitasatospora_viridis_DSM_44826.fa', 'Streptomyces_sp_1.fna', 'Other.fna']
# Only the first two genomes passed filtering
BLASTOUT = pd.DataFrame({'sseqid': pd.Categorical(['0', '1', '1'])})


@pytest.mark.parametrize(
    'given,expected',
    [
        ('', ''),
        ('Kitasatospora_viridis_DSM_44826', 'Kitasatospora_viridis_DSM_44826'),
        ('Kitasatospora_viridis_DSM_44826.fa', 'Kitasatospora_viridis_DSM_44826'),
        (
            '/Volumes/Seagate/Genomes/Kitasatospora_viridis_DSM_44826.fa',
            'Kitasatospora_viridis_DSM_44826',
        ),
        (
            'Kitasatospora_viridis_DSM_44826.fa, Streptomyces_sp_1',
            'Kitasatospora_viridis_DSM_44826,Streptomyces_sp_1',
        ),
    ],
)
def test_check_outgroup(given, expected):
    assert check_outgroup(given, LABELS, BLASTOUT) == expected


@pytest.mark.parametrize('given', ['Missing_genome', 'Other.fna'])
def test_check_outgroup_not_in_analysis(given):
    with pytest.raises(SystemExit) as err:
        check_outgroup(given, LABELS, BLASTOUT)
    assert err.value.code == 78
