import fidimag


def test_citation():
    s = fidimag.citation()
    assert type(s) == str
    assert 'fidimag' in s.lower()

    sb = fidimag.citation(bibtex=True)
    assert sb[0] == '@'
    assert len(sb.split()) > 1
    assert 'fidimag' in sb.lower()
           
