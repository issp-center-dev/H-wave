"""Regression guard for issue #147: in the Japanese manual, an inline role or
literal whose closing backtick abuts a full-width (East-Asian wide/fullwidth)
character is not recognised by docutils, so the raw markup is rendered instead
of the formula -- and one unclosed role throws off backtick pairing for the
rest of the paragraph. Sphinx does not warn about it (the build succeeds), so a
source-level check is the only cheap guard.

The fix is the RST escaped space ``\\ `` after the closing backtick(s); this
test fails if any un-escaped abutment reappears in ``docs/ja``.
"""
import os
import re
import unicodedata
import unittest

_HERE = os.path.dirname(os.path.abspath(__file__))
_JA_ROOT = os.path.normpath(os.path.join(_HERE, "..", "docs", "ja"))

#: an inline role ``:name:`...` `` or an inline literal ``` ``...`` ``` (no
#: embedded backtick); the group captures up to and including the closing
#: backtick(s) so the following character can be inspected
_SPAN = re.compile(r":[A-Za-z_+.\-]+:`[^`\n]+`|``[^`\n]+``")


def _is_wide(ch):
    return unicodedata.east_asian_width(ch) in ("W", "F")


def _offenders(text):
    """Line/column of every inline span whose close abuts a wide character."""
    hits = []
    for lineno, line in enumerate(text.split("\n"), 1):
        for m in _SPAN.finditer(line):
            nxt = line[m.end():m.end() + 1]
            if nxt and _is_wide(nxt):
                hits.append((lineno, line.strip()[:100]))
    return hits


class TestJapaneseDocsInlineMarkup(unittest.TestCase):
    def test_no_inline_markup_abuts_a_fullwidth_character(self):
        if not os.path.isdir(_JA_ROOT):
            self.skipTest("docs/ja not present")
        problems = []
        for dp, _, fns in os.walk(_JA_ROOT):
            for fn in sorted(fns):
                if not fn.endswith(".rst"):
                    continue
                p = os.path.join(dp, fn)
                with open(p, encoding="utf-8") as fh:
                    hits = _offenders(fh.read())
                for lineno, snippet in hits:
                    problems.append("{}:{}: {}".format(
                        os.path.relpath(p, _JA_ROOT), lineno, snippet))
        self.assertEqual(
            problems, [],
            "inline markup abuts a full-width character (issue #147); insert an "
            "RST escaped space '\\ ' after the closing backtick(s):\n"
            + "\n".join(problems))


if __name__ == "__main__":
    unittest.main()
