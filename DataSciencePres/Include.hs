#!/usr/bin/env runhaskell
import Text.Pandoc.JSON

doInclude :: Block -> IO Block
doInclude cb@(CodeBlock ("verbatim", classes, namevals) contents) =
  case lookup "include" namevals of
       Just f     -> return . (\x -> Para [Math DisplayMath x]) =<< readFile f
       Nothing    -> return cb
doInclude cb@(CodeBlock (id, classes, namevals) contents) =
  case lookup "include" namevals of
       Just f     -> return . (CodeBlock (id, classes, namevals)) =<< readFile f
       Nothing    -> return cb
doInclude x = return x

main :: IO ()
main = toJSONFilter doInclude
