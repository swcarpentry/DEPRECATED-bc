# OOP Exercises
Exercises by Orion Buske with thanks to Jon Pipitone.

--- 

Orion is excited at the prospect of having sushi for lunch tomorrow, so this
seems like a perfect opportunity to practice object-oriented programming. My
apologies in advance for the abundant over-simplifications and high likelihood
of other mistakes.

Let's start with a little bit of background. Sushi is, to quote Wikipedia, "a
Japanese dish consisting of cooked vinegared rice which is commonly topped with
other ingredients, such as fish or other seafood, or put into rolls." The two
most popular forms of sushi, with a few sub-types, are:

1. **Nigiri**- a hand-formed ball of rice topped with something tasty
  * Gunkanmaki - with a loose or soft topping that is held in place with a
    strip of seaweed
  * Temarizushi - where the topping is just pressed into a ball of rice

2. **Maki**- one or more tasty things rolled up in seaweed and rice
  * Futomaki - seaweed on the outside, usually vegetarian
  * Temaki - cone-shaped seaweed filled with rice and tasty things
  * Uramaki - rice on the outside

Luckily, this sort of hierarchical structure lends itself nicely to Classes and
inheritance. If you have been in a sushi restaurant before, you know how often
there are typos in the English descriptions. We are going to write a simple
program that a sushi restaurant owner could (theoretically) use to create a
menu, complete with English translation and Japanese transliteration (but not
actual Japanese, forgive me). 

## 1

To start, create an `Ingredient` class that inherits from `object`. The
constructor should accept two strings as arguments, `japanese` and `english`,
that correspond to the Japanese transliteration and English translation. The
`english` argument should be optional, and should default to the value of
`japanese` if not supplied (just like on menus, where some ingredients aren't
translated and you're left to wonder hopelessly). The value of both should be
saved as members of the Ingredient class.

**Hint:**
> ```python
> def __init__(self, japanese, english=None):
>     assert japanese
>     if english is None:
>         english = japanese
>     ...
> ```


## 2

Add to the `Ingredient` class two methods: `__str__(self)` and
`to_english(self)`. Both methods must return a string, and the `__str__` method
is what gets called when you print an object or "cast" it to a string. We will
have `__str__` return the Japanese name of the ingredient, and `to_english` will
return the English translation.

**Hint:**
> ```python
> def __str__(self):
>     return self.japanese
> 
> ```

I have created (with data from http://www.bento.com/sushivoc.html) a file of
common sushi ingredients, [oop_sushi.txt](oop_sushi.txt). The first column is the
transliteration, the second column is the translation, if available (I
selectively removed a few). 



## 3

Write a function, `read_ingredients`, that accepts an opened file object as its
argument and returns a list of `Ingredient` objects.

Try calling this function on the sushi_terms.txt in an `if __name__ ==
'__main__'` block and printing the first few ingredients to make sure it works.


## 4

Now, create a `Sushi` class that inherits from `object`. `Sushi`s It should have
a constructor that accepts a list of `Ingredient` objects.

## 5


Next, add a `__str__(self)` method. This method must return a string.  The
string should contain the Japanese representation of all the ingredients, but
the string itself should be in proper English so, for example, "buri", "buri and
tsubugai", and "buri, tsubugai, and kanpachi" are the correct way to print one,
two, or three ingredients, respectively. Do not just join the ingredients with
commas.

**Hint:**

> Since all the ingredients are `Ingredient` objects, you can just turn them
> into strings to get their Japanese representation.


**Hint:**

> There are three cases: 1, 2, or 3+ items. It's okay to handle them separately.


## 6

Next, add a loop to your `__main__` block that prompts the user for a menu item
and reads a line from `sys.stdin`. Provide a command for the user to quit (and
tell them what it is). For now, expect the user to just type one or more
ingredients on a line. You can use the built-in function `raw_input()` for this.

You should then parse the ingredients, find the appropriate `Ingredient`
objects, create a `Sushi` object, and print it in response. For example:

    
    Enter your sushi ('QUIT' to exit): <strong>unagi fugu ika sake</strong>
    unagi, fugu, ika, and sake


You may need to review
[dictionaries](http://software-carpentry.org/4_0/setdict/dict/) for this
exercise.



## 7

Now, add another method to the Sushi class, `to_english(self)`, which should
return the English translation for the Sushi object. Thus, it should return a
similar string as the `__str__` method, but with English ingredients instead of
Japanese ones. Do not call `__str__` and translate its string. Since you were
given the ingredients initially, just use their `to_english` methods, format
them correctly with commas and "and"s, and return that. Since both `to_english`
and `__str__` have to format their ingredients in the same way, you might want
to create a helper method that formats a list of ingredients (regardless of
their language).

You should now also print the result of calling `to_english` on the `Sushi`
objects you make at the user's request. Thus:
    
    Enter your sushi ('QUIT' to exit): <strong>unagi fugu ika sake</strong>
    unagi, fugu, ika, and sake
    eel, fugu, squid, and salmon


## 8

Now let's add a `Maki` class that inherits from `Sushi`. Everything will be the
same, except instead of just printing the ingredients, we want to print
something more descriptive. Let's have its `__str__` and `to_english` methods
return a string of the form: `[ingredients] rolled in [rice] and [seaweed]`,
where `[ingredients]` is our grammatical list of ingredients, and `[rice]` and
`[seaweed]` are two other ingredients that will be consistent across all sushi
types, but you should be sure to use the correct language at the correct time,
like other ingredients. However, these ingredients won't be specified in the
list of ingredients; they are implied by the type of sushi! You can create
constants for these ingredients or handle them in some other way. I did the
following:

```python    
    RICE = Ingredient('su-meshi', 'sushi rice')
    SEAWEED = Ingredient('nori', 'seaweed')
```


## 9

Now, revise the `__main__` block so that if someone enters "unagi fugu" or
"unagi fugu sushi" we consider it to be general sushi and create an appropriate
`Sushi` object. However, if the last word was "maki", we should create a `Maki`
object instead. You should do this in a way that is very easy to extend, because
there are going to be many more of these. As a general rule, we'll expect the
user to enter a number of Japanese ingredients, possibly following by a sushi
type. If no sushi type is specified, we should default to the base class,
otherwise we should use the type the user specified.

**Hint:**

```pythohn
    types = {'sushi': Sushi,
             'maki': Maki}
    if words[-1] not in types:
        Sushi(words)
    else:
        ...
```


## 10

Wonderful! We have a few more kinds of sushi to add, though. Futomaki, Temaki,
and Uramaki are all types of Maki, and all should inherit from it. Their
respective format strings should be of the following sort:

- **Futomaki:** "[ingredients] rolled in [rice] and [seaweed], with [seaweed]
  facing out"
- **Temaki:**   "cone of [seaweed] filled with [rice] and [ingredients]"
- **Uramaki:**  "[ingredients] rolled in [seaweed] and [rice], with [rice]
  facing out"

You may find the notion of a format string useful in this endeavor. For
instance, if you have the following string:

```python
    my_str = "%(ingredients)s rolled in %(rice)s and %(seaweed)s, with %(seaweed)s facing out"
```

Then you can do the following:

```python    
    >>> ingredient_str = "yummy things"
    >>> print my_str % {'rice': RICE, 'seaweed': SEAWEED, 'ingredients': ingredient_str}
    yummy things rolled in su-meshi and nori, with nori facing out
``

This is in fact quite powerful, because you can include rice and seaweed even if
they don't occur in the format string! Given this knowledge, you should try to
rewrite the `Sushi` base class so that it formats a member variable with a
dictionary of 'rice', 'seaweed', and 'ingredients'. Then, any child class need
only change their value of this member and everything works. For example:


```python
    class Futomaki(Maki):
        description = "%(ingredients)s rolled in %(rice)s and %(seaweed)s, with %(seaweed)s facing out"
    
    class Temaki(Maki):
        description = "cone of %(seaweed)s filled with %(rice)s and %(ingredients)s"
```

Make sure this works for both Japanese and English strings, and make sure you've
added these new Maki types to your `__main__` block so that the following work:
    
    Enter your sushi ('QUIT' to exit): <strong>unagi ohyo uramaki</strong>
    unagi and ohyo rolled in nori and su-meshi, with su-meshi facing out
    eel and halibut rolled in seaweed and sushi rice, with sushi rice facing out
    
    Enter your sushi ('QUIT' to exit): <strong>ikura temaki</strong>
    cone of nori filled with su-meshi and ikura
    cone of seaweed filled with sushi rice and salmon roe



## 11

Almost done. One last set of sushi classes to add. Add a Nigiri class that
inherits from Sushi, and Gunkanmaki and Temarizushi classes that inherit from
Nigiri. Since Nigiri usually only has one topping, you should take advantage of
inheritance to make sure this is true for all such sushi by checking that this
is the case in Nigiri's `__init__` method. If you run into an error, raise an
InvalidSushiError (you will have to define one; Python's libraries aren't quite
that complete). Don't forget to call it's parent's init method as well. Their
descriptions are as follows:
    
- **Nigiri:** "hand-formed [rice] topped with [ingredients]"
- **Gunkanmaki:** "[ingredient] on [rice] wrapped in a strip of [seaweed]"
- **Temarizushi:** "[ingredients] pressed into a ball of [rice]"

**Hint:**

```python
    class InvalidSushiError(Exception):
        pass
```

As a final test, the following example should work:

    
    Enter your sushi ('QUIT' to exit): <strong>fugu ohyo ika unagi</strong>
    fugu, ohyo, ika, and unagi
    fugu, halibut, squid, and eel
    
    Enter your sushi ('QUIT' to exit): <strong>fugu ohyo ika unagi sushi</strong>
    fugu, ohyo, ika, and unagi
    fugu, halibut, squid, and eel
    
    Enter your sushi ('QUIT' to exit): <strong>ika sake gunkanmaki</strong>
    Traceback (most recent call last):
        ...
    __main__.InvalidSushiError: Nigiri has only one topping



