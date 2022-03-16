import picmistandard
import unittest
import typing


class Test_ClassWithInit(unittest.TestCase):
    class MockClass(picmistandard.base._ClassWithInit):
        # note: refer to .base b/c class name with _ will not be exposed
        mandatory_attr: typing.Any
        name = ""
        optional = None
        _protected = 1

    class MockCheckTracer(picmistandard.base._ClassWithInit):
        """
        used to demonstrate the check interface
        """
        check_pass: bool = True
        check_counter = 0
        must_be_str: str = ""

        # used as static constant (though they dont actually exist in python)
        ERRORMSG = "apples-hammer-red"

        def _check(self) -> None:
            self.check_counter += 1
            # note: assign a specific message to assert for this exact
            # exception in tests
            assert self.check_pass, self.ERRORMSG

    def setUp(self):
        picmistandard.register_codename("mockpic")

    def test_arguments_used(self):
        """init sets provided args to attrs"""
        d = self.MockClass(mandatory_attr=None,
                                  name="n",
                                  optional=17)
        self.assertEqual(None, d.mandatory_attr)
        self.assertEqual("n", d.name)
        self.assertEqual(17, d.optional)

    def test_defaults(self):
        """if not given, defaults are used"""
        d = self.MockClass(mandatory_attr=42)
        self.assertEqual("", d.name)
        self.assertEqual(None, d.optional)

    def test_unkown_rejected(self):
        """unknown names are rejected"""
        with self.assertRaisesRegex(NameError, ".*blabla.*"):
            self.MockClass(mandatory_attr=1,
                                  blabla="foo")

    def test_codespecific(self):
        """arbitrary attrs for code-specific args used"""
        # args beginning with mockpic_ must be accepted
        d1 = self.MockClass(mandatory_attr=2,
                                   mockpic_foo="bar",
                                   mockpic_baz="xyzzy",
                                   mockpic=1,
                                   mockpic_=3)
        self.assertEqual("bar", d1.mockpic_foo)
        self.assertEqual("xyzzy", d1.mockpic_baz)
        self.assertEqual(1, d1.mockpic)
        self.assertEqual(3, d1.mockpic_)

        # _ separator is required:
        with self.assertRaisesRegex(NameError, ".*mockpicno_.*"):
            self.MockClass(mandatory_attr=2,
                                  mockpicno_="None")

        # args from other supported codes are still accepted
        d2 = self.MockClass(mandatory_attr=None,
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
            self.MockClass()

        # ok:
        d = self.MockClass(mandatory_attr="x")
        self.assertEqual("x", d.mandatory_attr)

    def test_typechecks(self):
        """typechecks only in check()"""
        class WithTypecheck(picmistandard.base._ClassWithInit):
            attr: str
            num: int = 0

        w = WithTypecheck(attr="d", num=2)

        # can overwrite vars, but then check fails
        w.attr = None
        with self.assertRaisesRegex(TypeError, ".*str.*"):
            w.check()

        # also checks in constructor:
        with self.assertRaisesRegex(TypeError, ".*str.*"):
            WithTypecheck(attr=7283)

        with self.assertRaisesRegex(TypeError, ".*int.*"):
            WithTypecheck(attr="", num=123.3)

    def test_protected(self):
        """protected args may *never* be accessed"""
        with self.assertRaisesRegex(NameError, ".*_protected.*"):
            self.MockClass(mandatory_attr=1,
                                  _protected=42)

        # though, *technically speaking*, it can be assigned
        d = self.MockClass(mandatory_attr=1)
        # ... this is evil, never do this!
        d._protected = 3
        self.assertEqual(3, d._protected)

    def test_check_basic(self):
        """simple demonstration of check() interface"""
        # passes
        check_tracer = self.MockCheckTracer()
        check_tracer.check()

        # make check() fail:
        check_tracer.check_pass = False
        with self.assertRaisesRegex(AssertionError,
                                    self.MockCheckTracer.ERRORMSG):
            check_tracer.check()

        with self.assertRaisesRegex(AssertionError,
                                    self.MockCheckTracer.ERRORMSG):
            self.MockCheckTracer(check_pass=False)

    def test_empty(self):
        """empty object works"""
        class MockEmpty(picmistandard.base._ClassWithInit):
            pass

        # both just pass
        empty = MockEmpty()
        empty.check()

    def test_check_optional(self):
        """implementing check() is not required"""
        class MockNoCheck(picmistandard.base._ClassWithInit):
            attr = 3

        no_check = MockNoCheck()
        # method exists & passes -- no matter the attribute value
        for value in [1, None, {}, [], ""]:
            no_check.attr = value
            no_check.check()

    def test_check_in_init(self):
        """check called from constructor"""
        check_tracer = self.MockCheckTracer()
        # counter is already one
        self.assertEqual(1, check_tracer.check_counter)

        # ... even if its default is zero
        self.assertEqual(0, check_tracer.__class__.__dict__["check_counter"])

        # one more call -> counter increased by one
        check_tracer.check()
        self.assertEqual(2, check_tracer.check_counter)

    def test_default_invalid_type(self):
        """raises if default variable has invalid type"""
        class MockInvalidDefaultType(picmistandard.base._ClassWithInit):
            my_str_attr: str = None

        with self.assertRaisesRegex(TypeError, ".*default.*my_str_attr.*"):
            MockInvalidDefaultType()

    def test_check_order(self):
        """_check() is only called if typechecks pass"""
        check_tracer = self.MockCheckTracer()

        cnt_old = check_tracer.check_counter

        # check will now fail *every time* when called
        check_tracer.check_pass = False

        # make type check break
        check_tracer.must_be_str = None
        with self.assertRaises(TypeError):
            check_tracer.check()

        # typecheck failed before _check() could be called
        # -> counter at old state
        self.assertEqual(cnt_old, check_tracer.check_counter)

        # when the type checks pass, _check is called (which fails)
        check_tracer.must_be_str = ""
        with self.assertRaisesRegex(AssertionError,
                                    self.MockCheckTracer.ERRORMSG):
            check_tracer.check()

        # counter increased
        self.assertEqual(cnt_old + 1, check_tracer.check_counter)

    def test_attribute_optional(self):
        """attributes can be (explicitly) made optional"""
        class MockOptionalAttrs(picmistandard.base._ClassWithInit):
            mandatory: str
            num_with_default: float = 3
            optional_name: typing.Optional[str] = None

        poa = MockOptionalAttrs(mandatory="", optional_name="foo")
        # optional_name can be set to none, and still passes:
        poa.optional_name = None
        poa.check()

        # but removing the mandatory arg raises:
        poa.mandatory = None
        with self.assertRaises(TypeError):
            # note: type error b/c NoneType != str
            poa.check()

    def test_dict_cast(self):
        """object can be cast to dict"""
        # this checks if a quirk of python is circumvented:
        # if you provide a default value for a class attribute, e.g.
        # >>> class MyClass:
        # ...   attr = 2
        # will result in
        # >>> 2 == MyClass().attr
        # True
        #
        # However, when casting the *instance* to a dictionary,
        # this variable will be missing:
        #
        # >>> myobject = MyClass()
        # >>> myobject.__dict__
        # {}
        #
        # only after explicitly setting the attribute will it become available:
        # >>> myobject.attr = myobject.attr
        # >>> myobject.__dict__
        # {'attr': 2}
        #
        # To aid serialization, we expect the __dict__ cast after construction
        # to be complete, e.g. for making class variables available in free
        # expressions.

        class MockWithDefaults(picmistandard.base._ClassWithInit):
            is_one = 1
            empty_str = ""
            empty_list = []

        mock_object = MockWithDefaults()
        self.assertEqual(mock_object.__dict__,
                         {"is_one": 1,
                          "empty_str": "",
                          "empty_list": []})
