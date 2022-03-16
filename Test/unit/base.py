import picmistandard
import unittest
import typing


class Test_ClassWithInit(unittest.TestCase):
    class PlaceholderClass(picmistandard.base._ClassWithInit):
        # note: refer to .base b/c class name with _ will not be exposed
        mandatory_attr: typing.Any
        name = ""
        optional = None
        _protected = 1

    def setUp(self):
        picmistandard.register_codename("placeholderpic")

    def test_arguments_used(self):
        """init sets provided args to attrs"""
        d = self.PlaceholderClass(mandatory_attr=None,
                                  name="n",
                                  optional=17)
        self.assertEqual(None, d.mandatory_attr)
        self.assertEqual("n", d.name)
        self.assertEqual(17, d.optional)

    def test_defaults(self):
        """if not given, defaults are used"""
        d = self.PlaceholderClass(mandatory_attr=42)
        self.assertEqual("", d.name)
        self.assertEqual(None, d.optional)

    def test_unkown_rejected(self):
        """unknown names are rejected"""
        with self.assertRaisesRegex(NameError, ".*blabla.*"):
            self.PlaceholderClass(mandatory_attr=1,
                                  blabla="foo")

    def test_codespecific(self):
        """arbitrary attrs for code-specific args used"""
        # args beginning with placeholderpic_ must be accepted
        d1 = self.PlaceholderClass(mandatory_attr=2,
                                   placeholderpic_foo="bar",
                                   placeholderpic_baz="xyzzy",
                                   placeholderpic=1,
                                   placeholderpic_=3)
        self.assertEqual("bar", d1.placeholderpic_foo)
        self.assertEqual("xyzzy", d1.placeholderpic_baz)
        self.assertEqual(1, d1.placeholderpic)
        self.assertEqual(3, d1.placeholderpic_)

        # _ separator is required:
        with self.assertRaisesRegex(NameError, ".*placeholderpicno_.*"):
            self.PlaceholderClass(mandatory_attr=2,
                                  placeholderpicno_="None")

        # args from other supported codes are still accepted
        d2 = self.PlaceholderClass(mandatory_attr=None,
                                   warpx_anyvar=1,
                                   warpx=2,
                                   warpx_=3,
                                   fbpic=4)
        self.assertEqual(None, d2.mandatory_attr)
        self.assertEqual(1, d2.warpx_anyvar)
        self.assertEqual(2, d2.warpx)
        self.assertEqual(3, d2.warpx_)
        self.assertEqual(4, d2.fbpic)

    def test_mandatory_enforced(self):
        """mandatory args must be given"""
        with self.assertRaisesRegex(RuntimeError, ".*mandatory_attr.*"):
            self.PlaceholderClass()

        # ok:
        d = self.PlaceholderClass(mandatory_attr="x")
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

    def test_protected(self):
        """protected args may *never* be accessed"""
        with self.assertRaisesRegex(NameError, ".*_protected.*"):
            self.PlaceholderClass(mandatory_attr=1,
                                  _protected=42)

        # though, *technically speaking*, it can be assigned
        d = self.PlaceholderClass(mandatory_attr=1)
        # ... this is evil, never do this!
        d._protected = 3
        self.assertEqual(3, d._protected)
