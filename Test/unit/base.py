import picmistandard
import unittest
import typing


class Test_ClassWithInit(unittest.TestCase):
    class DummyClass(picmistandard.base._ClassWithInit):
        # note: refer to .base b/c class name with _ will not be exposed
        mandatory_attr: typing.Any
        name = ""
        optional = None

    def setUp(self):
        picmistandard.register_codename("dummypic")

    def test_arguments_used(self):
        """init sets provided args to attrs"""
        d = DummyClass(mandatory_attr=None,
                       name="n",
                       optional=17)
        self.assertEqual(None, d.mandatory_attr)
        self.assertEqual("n", d.name)
        self.assertEqual(17, d.optional)

    def test_defaults(self):
        """if not given, defaults are used"""
        d = DummyClass(mandatory_attr=42)
        self.assertEqual("", d.name)
        self.assertEqual(None, d.optional)

    def test_unkown_rejected(self):
        """unknown names are rejected"""
        with self.assertRaisesRegex(NameError, ".*blabla.*"):
            DummyClass(mandatory_attr=1,
                       blabla="foo")

    def test_codespecific(self):
        """arbitrary attrs for code-specific args used"""
        # args beginning with dummypic_ must be accepted
        d1 = DummyClass(mandatory_attr=2,
                        dummypic_foo="bar",
                        dummypic_baz="xyzzy")
        self.assertEqual("bar", d1.dummypic_foo)
        self.assertEqual("xyzzy", d1.dummypic_baz)

        # _ separator is required:
        with self.assertRaisesRegex(NameError, ".*dummypicno_.*"):
            DummyClass(mandatory_attr=2,
                       dummypicno_="None")

        # args from other supported codes are still accepted
        d2 = DummyClass(mandatory_attr=None,
                        warpx_anyvar=1)
        self.assertEqual(None, d2.mandatory_attr)
        self.assertEqual(1, d2.warpx_anyvar)

    def test_mandatory_enforced(self):
        """mandatory args must be given"""
        with self.assertRaisesRegex(RuntimeError, ".*mandatory_attr.*"):
            DummyClass()

        # ok:
        d = DummyClass(mandatory_attr="x")
        self.assertEqual("x", d.mandatory_attr)

    def test_no_typechecks(self):
        """no typechecks, explicit type annotations are rejected"""
        class WithTypecheck(picmistandard.base._ClassWithInit):
            attr: str
            num: int = 0

        with self.assertRaises(SyntaxError):
            # must complain purely b/c typecheck is *there*
            # (even if it would enforceable)
            WithTypecheck(attr="d", num=2)
